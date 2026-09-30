import argparse
import csv
import json
import os
import pathlib
import platform
import subprocess
import sys
import time

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from sklearn.metrics import average_precision_score, roc_auc_score

PROJECT_ROOT = pathlib.Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

# Local NPZ loader replaces unrelated, unused gensim/Word2Vec dependency.
class DataLoader:
    @staticmethod
    def load_embeddings(path):
        data = np.load(path, allow_pickle=False)
        return data['ID'], data['embeddings']



def l2_normalize_rows(x: np.ndarray) -> np.ndarray:
    norms = np.linalg.norm(x, axis=1, keepdims=True)
    norms = np.where(norms == 0.0, 1.0, norms)
    return x / norms


def transform_weights(w: np.ndarray, method: str) -> np.ndarray:
    m = method.lower()
    if m == "none":
        return w
    if m == "log1p":
        return np.log1p(np.maximum(w, 0.0))
    if m == "sqrt":
        return np.sqrt(np.maximum(w, 0.0))
    raise ValueError(f"Unknown edge_weight_transform: {method}")


def build_edges(
    edge_file: str,
    weight_column: str,
    id_to_idx: dict[str, int],
    edge_weight_transform: str,
) -> tuple[np.ndarray, np.ndarray]:
    df = pd.read_csv(edge_file, delimiter="\t")
    required = {"GS_A_ID", "GS_B_ID", weight_column}
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError(f"Missing columns in edge_file: {missing}")

    rows = []
    for _, row in df.iterrows():
        a = str(row["GS_A_ID"])
        b = str(row["GS_B_ID"])
        if a not in id_to_idx or b not in id_to_idx:
            continue
        u = id_to_idx[a]
        v = id_to_idx[b]
        if u == v:
            continue
        if u > v:
            u, v = v, u
        rows.append((u, v, float(row[weight_column])))

    if not rows:
        raise ValueError("No aligned edges")

    edge_df = (
        pd.DataFrame(rows, columns=["u", "v", "w"])
        .groupby(["u", "v"], as_index=False)["w"]
        .max()
    )
    edge_df["w"] = transform_weights(edge_df["w"].to_numpy(), edge_weight_transform)
    edge_index = edge_df[["u", "v"]].to_numpy(dtype=np.int64)
    edge_weight = edge_df["w"].to_numpy(dtype=np.float32)
    return edge_index, edge_weight


def encode_pair_key(u: np.ndarray, v: np.ndarray, n: int) -> np.ndarray:
    return u.astype(np.int64) * np.int64(n) + v.astype(np.int64)


def build_normalized_adjacency(
    num_nodes: int,
    edge_index: np.ndarray,
    edge_weight: np.ndarray,
    device: torch.device,
) -> torch.Tensor:
    u = edge_index[:, 0]
    v = edge_index[:, 1]
    row = np.concatenate([u, v, np.arange(num_nodes, dtype=np.int64)])
    col = np.concatenate([v, u, np.arange(num_nodes, dtype=np.int64)])
    val = np.concatenate(
        [edge_weight, edge_weight, np.ones(num_nodes, dtype=np.float32)]
    )

    row_t = torch.from_numpy(row).to(device)
    col_t = torch.from_numpy(col).to(device)
    val_t = torch.from_numpy(val).to(device)

    deg = torch.zeros(num_nodes, dtype=torch.float32, device=device)
    deg.index_add_(0, row_t, val_t)
    deg_inv_sqrt = torch.pow(deg.clamp(min=1e-12), -0.5)
    norm_val = deg_inv_sqrt[row_t] * val_t * deg_inv_sqrt[col_t]

    idx = torch.stack([row_t, col_t], dim=0)
    return torch.sparse_coo_tensor(
        idx, norm_val, (num_nodes, num_nodes), device=device
    ).coalesce()


class GraphSageEncoder(nn.Module):
    def __init__(
        self, input_dim: int, hidden_dim: int, num_layers: int, dropout: float
    ) -> None:
        super().__init__()
        if num_layers < 1:
            raise ValueError("num_layers must be >=1")
        self.layers = nn.ModuleList()
        self.dropout = dropout
        d = input_dim
        for _ in range(num_layers - 1):
            self.layers.append(nn.Linear(d * 2, hidden_dim))
            d = hidden_dim
        self.out = nn.Linear(d * 2, hidden_dim)

    def forward(self, x: torch.Tensor, adj: torch.Tensor) -> torch.Tensor:
        h = x
        for layer in self.layers:
            neigh = torch.sparse.mm(adj, h)
            h = torch.cat([h, neigh], dim=1)
            h = F.relu(layer(h))
            h = F.dropout(h, p=self.dropout, training=self.training)
        neigh = torch.sparse.mm(adj, h)
        h = torch.cat([h, neigh], dim=1)
        return self.out(h)


def dot_decoder(z: torch.Tensor, edges: torch.Tensor) -> torch.Tensor:
    return (z[edges[:, 0]] * z[edges[:, 1]]).sum(dim=1)


def sample_negative_edges(num_nodes, n_samples, edge_key_set, rng):
    # Enumerating this small graph avoids duplicate negatives and infinite retries.
    u, v = np.triu_indices(num_nodes, k=1)
    keys = encode_pair_key(u, v, num_nodes)
    available = np.flatnonzero(np.array([int(k) not in edge_key_set for k in keys]))
    if n_samples > len(available):
        raise ValueError('Insufficient distinct non-edges for requested split')
    picked = rng.choice(available, size=n_samples, replace=False)
    return np.column_stack([u[picked], v[picked]])


def evaluate_auc_ap(
    model: GraphSageEncoder,
    x: torch.Tensor,
    adj: torch.Tensor,
    pos_edges: np.ndarray,
    neg_edges: np.ndarray,
    device: torch.device,
) -> tuple[float, float]:
    model.eval()
    with torch.no_grad():
        z = model(x, adj)
        pos_t = torch.from_numpy(pos_edges).to(device)
        neg_t = torch.from_numpy(neg_edges).to(device)
        pos_scores = torch.sigmoid(dot_decoder(z, pos_t)).cpu().numpy()
        neg_scores = torch.sigmoid(dot_decoder(z, neg_t)).cpu().numpy()
    y_true = np.concatenate([np.ones_like(pos_scores), np.zeros_like(neg_scores)])
    y_score = np.concatenate([pos_scores, neg_scores])
    return float(roc_auc_score(y_true, y_score)), float(
        average_precision_score(y_true, y_score)
    )


def evaluate_loss_auc_ap(
    model: GraphSageEncoder,
    x: torch.Tensor,
    adj: torch.Tensor,
    pos_edges: np.ndarray,
    neg_edges: np.ndarray,
    device: torch.device,
) -> tuple[float, float, float]:
    model.eval()
    with torch.no_grad():
        z = model(x, adj)
        pos_t = torch.from_numpy(pos_edges).to(device)
        neg_t = torch.from_numpy(neg_edges).to(device)
        pos_logits = dot_decoder(z, pos_t)
        neg_logits = dot_decoder(z, neg_t)
        loss = F.binary_cross_entropy_with_logits(pos_logits, torch.ones_like(pos_logits))
        loss = loss + F.binary_cross_entropy_with_logits(neg_logits, torch.zeros_like(neg_logits))

        pos_scores = torch.sigmoid(pos_logits).cpu().numpy()
        neg_scores = torch.sigmoid(neg_logits).cpu().numpy()

    y_true = np.concatenate([np.ones_like(pos_scores), np.zeros_like(neg_scores)])
    y_score = np.concatenate([pos_scores, neg_scores])
    auc = float(roc_auc_score(y_true, y_score))
    ap = float(average_precision_score(y_true, y_score))
    return float(loss.item()), auc, ap


def count_trainable_params(model: nn.Module) -> int:
    return int(sum(p.numel() for p in model.parameters() if p.requires_grad))


def maybe_git_sha(project_root: pathlib.Path) -> str:
    try:
        proc = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=str(project_root),
            check=False,
            capture_output=True,
            text=True,
        )
    except Exception:
        return ""
    if proc.returncode != 0:
        return ""
    return proc.stdout.strip()


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Holdout link prediction evaluation for GraphSAGE-style unsupervised training."
    )
    parser.add_argument("--edge_file", required=True)
    parser.add_argument("--feature_npz", required=True)
    parser.add_argument("--weight_column", default="NLOGPMF")
    parser.add_argument(
        "--edge_weight_transform", default="log1p", choices=["none", "log1p", "sqrt"]
    )
    parser.add_argument("--hidden_dim", type=int, default=256)
    parser.add_argument("--num_layers", type=int, default=2)
    parser.add_argument("--dropout", type=float, default=0.1)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--weight_decay", type=float, default=0.0)
    parser.add_argument("--epochs", type=int, default=120)
    parser.add_argument("--eval_every", type=int, default=10)
    parser.add_argument("--patience", type=int, default=0)
    parser.add_argument("--val_ratio", type=float, default=0.1)
    parser.add_argument("--test_ratio", type=float, default=0.1)
    parser.add_argument("--pos_samples_per_epoch", type=int, default=50000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--split_seed",
        type=int,
        default=-1,
        help="Seed for edge splitting + fixed val/test negative sampling. Defaults to --seed.",
    )
    parser.add_argument("--device", default="auto", choices=["auto", "cpu", "cuda"])
    parser.add_argument("--output_csv", required=True)
    parser.add_argument("--output_npz", default="")
    args = parser.parse_args()

    if args.split_seed < 0:
        args.split_seed = args.seed

    if args.val_ratio <= 0.0 or args.test_ratio <= 0.0:
        raise ValueError("val_ratio and test_ratio must be > 0")
    if args.val_ratio + args.test_ratio >= 1.0:
        raise ValueError("val_ratio + test_ratio must be < 1.0")

    np.random.seed(args.seed)
    torch.manual_seed(args.seed)
    rng_train = np.random.default_rng(args.seed)
    rng_split = np.random.default_rng(args.split_seed)

    if args.device == "auto":
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    else:
        device = torch.device(args.device)

    ids, emb = DataLoader.load_embeddings(args.feature_npz)
    ids = [str(x) for x in ids]
    x_np = l2_normalize_rows(np.asarray(emb, dtype=np.float32))
    num_nodes = x_np.shape[0]
    id_to_idx = {gid: i for i, gid in enumerate(ids)}

    edge_index, edge_weight = build_edges(
        args.edge_file, args.weight_column, id_to_idx, args.edge_weight_transform
    )
    n_edges = edge_index.shape[0]
    perm = rng_split.permutation(n_edges)
    n_test = max(1, int(args.test_ratio * n_edges))
    n_val = max(1, int(args.val_ratio * n_edges))
    if n_test + n_val >= n_edges:
        raise ValueError(
            f"Not enough edges for train split: n_edges={n_edges} n_test={n_test} n_val={n_val}"
        )
    test_idx = perm[:n_test]
    val_idx = perm[n_test : n_test + n_val]
    train_idx = perm[n_test + n_val :]

    train_edges = edge_index[train_idx]
    train_w = edge_weight[train_idx]
    val_edges = edge_index[val_idx]
    test_edges = edge_index[test_idx]

    edge_keys_all = set(
        encode_pair_key(edge_index[:, 0], edge_index[:, 1], num_nodes).tolist()
    )
    edge_keys_train = set(
        encode_pair_key(train_edges[:, 0], train_edges[:, 1], num_nodes).tolist()
    )
    val_neg = sample_negative_edges(num_nodes, len(val_edges), edge_keys_all, rng_split)
    val_neg_keys = set(encode_pair_key(val_neg[:, 0], val_neg[:, 1], num_nodes).tolist())
    test_neg = sample_negative_edges(num_nodes, len(test_edges), edge_keys_all | val_neg_keys, rng_split)
    test_neg_keys = set(encode_pair_key(test_neg[:, 0], test_neg[:, 1], num_nodes).tolist())
    train_negative_exclusions = edge_keys_all | val_neg_keys | test_neg_keys
    assert not (val_neg_keys & test_neg_keys)
    assert not ((val_neg_keys | test_neg_keys) & edge_keys_all)
    assert len(val_neg_keys) == len(val_edges) and len(test_neg_keys) == len(test_edges)
    pathlib.Path(args.output_csv).parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.output_csv + '.splits.npz', ID=np.asarray(ids, dtype=str),
                        train_positive=train_edges, val_positive=val_edges, test_positive=test_edges,
                        val_negative=val_neg, test_negative=test_neg)
    # Baseline uses the same held-out pairs and no graph training.
    yp = np.r_[np.ones(len(test_edges)), np.zeros(len(test_neg))]
    pairs = np.vstack([test_edges, test_neg])
    cosine = (x_np[pairs[:, 0]] * x_np[pairs[:, 1]]).sum(axis=1)
    baseline = {'semantic_cosine_test_auc': float(roc_auc_score(yp, cosine)),
                'semantic_cosine_test_ap': float(average_precision_score(yp, cosine)),
                'split_seed': args.split_seed, 'negative_splits_disjoint': True,
                'train_negative_excludes_all_true_edges_and_heldout_negatives': True}
    pathlib.Path(args.output_csv + '.baseline.json').write_text(json.dumps(baseline, indent=2))

    adj = build_normalized_adjacency(num_nodes, train_edges, train_w, device)
    x = torch.from_numpy(x_np).to(device)

    model = GraphSageEncoder(
        x_np.shape[1], args.hidden_dim, args.num_layers, args.dropout
    ).to(device)
    optimizer = torch.optim.Adam(
        model.parameters(),
        lr=args.lr,
        weight_decay=args.weight_decay,
    )

    out_dir = os.path.dirname(args.output_csv)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    print(
        f"device={device} nodes={num_nodes} edges(train/val/test)={len(train_edges)}/{len(val_edges)}/{len(test_edges)}",
        flush=True,
    )
    t0 = time.time()
    best_auc = -np.inf
    best_state: dict[str, torch.Tensor] | None = None
    best_epoch = 0
    best_metrics: dict[str, float] = {}
    no_improve = 0
    with open(args.output_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "epoch",
                "train_loss",
                "val_loss",
                "val_auc",
                "val_ap",
                "test_auc",
                "test_ap",
                "pos_samples",
                "elapsed_sec",
            ]
        )

        for epoch in range(1, args.epochs + 1):
            model.train()
            optimizer.zero_grad()
            z = model(x, adj)

            n_pos = min(args.pos_samples_per_epoch, len(train_edges))
            pos_idx = rng_train.choice(len(train_edges), size=n_pos, replace=False)
            pos = train_edges[pos_idx]
            neg = sample_negative_edges(num_nodes, n_pos, train_negative_exclusions, rng_train)

            pos_t = torch.from_numpy(pos).to(device)
            neg_t = torch.from_numpy(neg).to(device)

            pos_logits = dot_decoder(z, pos_t)
            neg_logits = dot_decoder(z, neg_t)

            loss = F.binary_cross_entropy_with_logits(
                pos_logits, torch.ones_like(pos_logits)
            )
            loss = loss + F.binary_cross_entropy_with_logits(
                neg_logits, torch.zeros_like(neg_logits)
            )
            loss.backward()
            optimizer.step()

            do_eval = (
                (epoch == 1) or (epoch % args.eval_every == 0) or (epoch == args.epochs)
            )
            if do_eval:
                val_loss, auc, ap = evaluate_loss_auc_ap(
                    model, x, adj, val_edges, val_neg, device
                )
                # Evaluate test only once, after validation-selected checkpoint.
                test_auc, test_ap = float('nan'), float('nan')
                elapsed = time.time() - t0
                w.writerow(
                    [
                        epoch,
                        float(loss.item()),
                        float(val_loss),
                        auc,
                        ap,
                        test_auc,
                        test_ap,
                        n_pos,
                        elapsed,
                    ]
                )
                f.flush()
                print(
                    f"epoch={epoch:03d} loss={loss.item():.4f} val_loss={val_loss:.4f} "
                    f"val_auc={auc:.4f} val_ap={ap:.4f} test_auc={test_auc:.4f} test_ap={test_ap:.4f}",
                    flush=True,
                )
                improved = auc > (best_auc + 1e-8)
                if improved:
                    best_auc = auc
                    best_epoch = epoch
                    best_state = {
                        k: v.detach().cpu().clone() for k, v in model.state_dict().items()
                    }
                    best_metrics = {
                        "train_loss": float(loss.item()),
                        "val_loss": float(val_loss),
                        "val_auc": float(auc),
                        "val_ap": float(ap),
                        "test_auc": float(test_auc),
                        "test_ap": float(test_ap),
                    }
                    no_improve = 0
                elif args.patience > 0:
                    no_improve += 1
                    if no_improve >= args.patience:
                        print(
                            f"early_stop at epoch={epoch} (patience={args.patience}, best_auc={best_auc:.4f})",
                            flush=True,
                        )
                        break

    # Final evaluation at best checkpoint (by val AUC).
    final_metrics: dict[str, float] = {}
    if best_state is not None:
        model.load_state_dict(best_state)
    final_val_loss, final_val_auc, final_val_ap = evaluate_loss_auc_ap(
        model, x, adj, val_edges, val_neg, device
    )
    final_test_loss, final_test_auc, final_test_ap = evaluate_loss_auc_ap(
        model, x, adj, test_edges, test_neg, device
    )
    final_metrics.update(
        {
            "val_loss": float(final_val_loss),
            "val_auc": float(final_val_auc),
            "val_ap": float(final_val_ap),
            "test_loss": float(final_test_loss),
            "test_auc": float(final_test_auc),
            "test_ap": float(final_test_ap),
        }
    )

    if args.output_npz:
        out_emb_dir = os.path.dirname(args.output_npz)
        if out_emb_dir:
            os.makedirs(out_emb_dir, exist_ok=True)
        model.eval()
        with torch.no_grad():
            z_final = model(x, adj).cpu().numpy().astype(np.float32)
        np.savez(args.output_npz, ID=np.asarray(ids, dtype=str), embeddings=z_final)
        print(f"saved embeddings: {args.output_npz}", flush=True)

    elapsed_total = time.time() - t0
    summary = {
        "script": str(pathlib.Path(__file__).resolve()),
        "git_sha": maybe_git_sha(PROJECT_ROOT),
        "env": {
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "torch": torch.__version__,
            "cuda_available": bool(torch.cuda.is_available()),
            "cuda_version": getattr(torch.version, "cuda", None),
            "device": str(device),
            "cuda_device_name": torch.cuda.get_device_name(0) if device.type == "cuda" else "",
        },
        "data": {
            "edge_file": args.edge_file,
            "feature_npz": args.feature_npz,
            "num_nodes": int(num_nodes),
            "input_dim": int(x_np.shape[1]),
            "n_edges_total": int(n_edges),
            "n_edges_train": int(len(train_edges)),
            "n_edges_val": int(len(val_edges)),
            "n_edges_test": int(len(test_edges)),
        },
        "config": {
            "weight_column": args.weight_column,
            "edge_weight_transform": args.edge_weight_transform,
            "hidden_dim": int(args.hidden_dim),
            "num_layers": int(args.num_layers),
            "dropout": float(args.dropout),
            "lr": float(args.lr),
            "weight_decay": float(args.weight_decay),
            "epochs": int(args.epochs),
            "eval_every": int(args.eval_every),
            "patience": int(args.patience),
            "val_ratio": float(args.val_ratio),
            "test_ratio": float(args.test_ratio),
            "pos_samples_per_epoch": int(args.pos_samples_per_epoch),
            "seed": int(args.seed),
            "split_seed": int(args.split_seed),
        },
        "model": {"n_params_trainable": count_trainable_params(model)},
        "best": {"epoch": int(best_epoch), **best_metrics},
        "final_at_best": final_metrics,
        "outputs": {
            "metrics_csv": args.output_csv,
            "summary_json": str(pathlib.Path(args.output_csv).with_suffix(".summary.json")),
            "embedding_npz": args.output_npz,
        },
        "runtime": {"elapsed_sec": float(elapsed_total)},
    }
    summary_path = pathlib.Path(args.output_csv).with_suffix(".summary.json")
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    print(f"wrote metrics: {args.output_csv}", flush=True)
    print(f"wrote summary: {summary_path}", flush=True)


if __name__ == "__main__":
    main()
