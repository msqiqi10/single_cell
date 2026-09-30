"""Create an auditable local patch of the inspected upstream GraphSAGE runner."""
from pathlib import Path
import json,hashlib,difflib
ROOT=Path(__file__).resolve().parents[1]
source=ROOT/'vendor/GOLDEN_GNN_gnn_linkpred_holdout.py'
old=source.read_text(); new=old
new=new.replace('from data_processing.data_loader import DataLoader', '''# Local NPZ loader replaces unrelated, unused gensim/Word2Vec dependency.
class DataLoader:
    @staticmethod
    def load_embeddings(path):
        data = np.load(path, allow_pickle=False)
        return data['ID'], data['embeddings']
''')
start=new.index('def sample_negative_edges(');end=new.index('\n\ndef evaluate_auc_ap',start)
new=new[:start]+'''def sample_negative_edges(num_nodes, n_samples, edge_key_set, rng):
    # Enumerating this small graph avoids duplicate negatives and infinite retries.
    u, v = np.triu_indices(num_nodes, k=1)
    keys = encode_pair_key(u, v, num_nodes)
    available = np.flatnonzero(np.array([int(k) not in edge_key_set for k in keys]))
    if n_samples > len(available):
        raise ValueError('Insufficient distinct non-edges for requested split')
    picked = rng.choice(available, size=n_samples, replace=False)
    return np.column_stack([u[picked], v[picked]])
''' + new[end:]
before='''    val_neg = sample_negative_edges(num_nodes, len(val_edges), edge_keys_all, rng_split)
    test_neg = sample_negative_edges(num_nodes, len(test_edges), edge_keys_all, rng_split)'''
after='''    val_neg = sample_negative_edges(num_nodes, len(val_edges), edge_keys_all, rng_split)
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
    pathlib.Path(args.output_csv + '.baseline.json').write_text(json.dumps(baseline, indent=2))'''
assert before in new;new=new.replace(before,after)
new=new.replace('neg = sample_negative_edges(num_nodes, n_pos, edge_keys_train, rng_train)',
                'neg = sample_negative_edges(num_nodes, n_pos, train_negative_exclusions, rng_train)')
before='''                _, test_auc, test_ap = evaluate_loss_auc_ap(
                    model, x, adj, test_edges, test_neg, device
                )'''
assert before in new;new=new.replace(before,'''                # Evaluate test only once, after validation-selected checkpoint.
                test_auc, test_ap = float('nan'), float('nan')''')
new=new.replace('ID=np.asarray(ids, dtype=object)', 'ID=np.asarray(ids, dtype=str)')
target=ROOT/'src/gnn_mouse_holdout.py';target.write_text(new)
patch=''.join(difflib.unified_diff(old.splitlines(True),new.splitlines(True),fromfile='upstream/gnn_linkpred_holdout.py',tofile='local/gnn_mouse_holdout.py'))
(ROOT/'provenance/gnn_adapter.patch').write_text(patch)
(ROOT/'provenance/gnn_adapter.json').write_text(json.dumps(dict(upstream_commit='a85654b6b853ddd60e59efa7f3bb769832d687b0',
    upstream_sha256=hashlib.sha256(old.encode()).hexdigest(),adapter_sha256=hashlib.sha256(new.encode()).hexdigest(),
    changes=['remove unused gensim dependency','unique and disjoint validation/test negative edges',
             'training negatives exclude ALL true edges including held-out positives and held-out negatives',
             'test evaluation only at final validation-selected checkpoint','export actual splits and semantic-only baseline',
             'non-pickled Unicode ID arrays'],architecture='upstream GraphSageEncoder unchanged; 256 hidden, 2 layers, dropout0.1, Adam1e-3,120 epochs'),indent=2))
compile(new,str(target),'exec');print('Adapter prepared and compiled',target)
