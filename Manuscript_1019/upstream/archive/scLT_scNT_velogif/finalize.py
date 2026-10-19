"""Validate delivered artifacts and produce a task-level report."""
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path('/task')


def main():
    validations = []
    statuses = []
    for dataset, basis, cluster, expected in [('scLT', 'cce', 'cell_type2', 3186),
                                                ('scNT', 'umap', 'cell_type', 1947)]:
        result = ROOT / 'results' / dataset
        summary = json.loads((result / 'run_summary.json').read_text())
        for item in summary['results']:
            statuses.append({'dataset': dataset, 'method': item['method'], 'status': item['status'],
                             'reason': item['message']})
        for method in ['topicvelo', 'bivi', 'deepkinet', 'cell2fate']:
            statuses.append({'dataset': dataset, 'method': method, 'status': 'excluded_missing_input',
                             'reason': 'Missing raw spliced/unspliced counts'})
        for filename in ['scvelo.sto.h5ad', 'scvelo.dyn.h5ad', 'sdevelo.h5ad']:
            path = result / filename
            if not path.exists():
                validations.append({'dataset': dataset, 'file': filename, 'status': 'missing'})
                continue
            obj = ad.read_h5ad(path)
            v = obj.layers['velocity']
            values = v.data if sparse.issparse(v) else np.asarray(v)
            assert obj.obs_names.is_unique
            assert obj.n_obs == expected
            assert cluster in obj.obs and 'X_' + basis in obj.obsm
            assert v.shape == obj.shape and np.isfinite(values).any()
            validations.append({'dataset': dataset, 'file': filename, 'status': 'validated',
                                'cells': obj.n_obs, 'genes': obj.n_vars,
                                'all_nan_velocity_genes': int(np.isnan(v.toarray() if sparse.issparse(v) else v).all(axis=0).sum()),
                                'nonfinite_stored_values': int((~np.isfinite(values)).sum()),
                                'bytes': path.stat().st_size})
        assert len(pd.read_csv(result / 'reference_evaluation/summary.csv')) == 3
        assert len(pd.read_csv(result / 'reference_evaluation/failures.csv')) == 0
        assert len(pd.read_csv(result / 'evals/Plot_failures.csv')) == 0
        assert len(list((result / 'figures').glob('*.svg'))) == 4
    pd.DataFrame(statuses).to_csv(ROOT / 'results/method_status.csv', index=False)
    (ROOT / 'records/final_validation.json').write_text(json.dumps(validations, indent=2))
    table = pd.read_csv(ROOT / 'results/reference_summary.csv')
    fig, ax = plt.subplots(figsize=(9, 5))
    table.pivot(index='method', columns='dataset', values='mean_cosine').plot.bar(ax=ax)
    ax.set_ylabel('Mean reference direction cosine')
    ax.set_xlabel('')
    ax.set_ylim(-1.05, 1.05)
    ax.axhline(0, color='grey', lw=.8)
    ax.tick_params(axis='x', rotation=15)
    ax.set_title('scLT and scNT (different reference systems)')
    fig.tight_layout()
    fig.savefig(ROOT / 'results/reference_summary.svg')
    fig.savefig(ROOT / 'results/reference_summary.png', dpi=150)
    plt.close(fig)
    lines = ['# 本次运行结果', '', '两套数据各完成 scVelo stochastic、dynamic 和 SDEvelo。', '',
             '| 数据集 | 方法 | 有效参考细胞 | 平均余弦 | 中位余弦 |',
             '|---|---|---:|---:|---:|']
    for _, r in table.iterrows():
        lines.append('| %s | %s | %d | %.5f | %.5f |' %
                     (r.dataset, r.method, r.valid_cells, r.mean_cosine, r.median_cosine))
    lines.extend(['', '参见 [跨数据集评价图](results/reference_summary.svg)、'
                  '[完整精度评价表](results/reference_summary.csv) 和 '
                  '[方法状态表](results/method_status.csv)。', '',
                  'TopicVelo、biVI、DeepKINET、cell2fate 因缺少原始 counts 未运行。',
                  '参考分数只统计坐标一致、有限且非零的速度，其他情况在逐细胞 CSV 中记录。',
                  'scLT 使用谱系参考的 CCE 方向，scNT 使用代谢参考的 UMAP 方向；'
                  '不同数据集的绝对分数不作为同一基准排名。',
                  '历史 CSV 保存在参考输入中，描述统计另存；本表全部来自本次新计算。'])
    lines.extend(['', '## 实际基因覆盖', '',
                  '| 数据集 | 文件 | 保留基因数 | velocity 全 NaN 基因数 |',
                  '|---|---|---:|---:|'])
    for row in validations:
        if row['status'] == 'validated':
            lines.append('| %s | %s | %d | %d |' % (row['dataset'], row['file'], row['genes'], row['all_nan_velocity_genes']))
    lines.extend(['', 'scVelo dynamic 对无法拟合的基因保留 NaN，投影时排除整列 NaN 基因。'
                  'SDEvelo 按其默认流程过滤基因并重新计算 moments，记录中的基因数反映实际输出。'])
    (ROOT / 'RESULTS.zh-CN.md').write_text('\n'.join(lines) + '\n')
    print(table.to_string(index=False))


if __name__ == '__main__':
    main()
