import os
import argparse
import anndata as ad
import numpy as np
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scvi
from sklearn.metrics import classification_report
from scarches.models.scpoli import scPoli
from scarches.dataset.trvae.data_handling import remove_sparsity

import warnings
warnings.filterwarnings('ignore')

class scPoliIntegration:
    def __init__(self, seed=42, cell_type_key="reanno", n_top_genes=2000, model_dir=None):
        """初始化scPoli整合分析类"""
        self.seed = seed
        self.cell_type_key = cell_type_key
        self.n_top_genes = n_top_genes
        self.model_dir = model_dir
        self.setup_environment()

    def setup_environment(self):
        """设置绘图和计算环境"""
        sc.settings.set_figure_params(dpi=100, frameon=False, figsize=(4, 4))
        plt.rcParams['figure.dpi'] = 100
        plt.rcParams['figure.figsize'] = (4, 4)
        sc.settings.seed = self.seed
        print(f"环境设置完成 - 细胞类型键: {self.cell_type_key}, 高变基因数: {self.n_top_genes}")

    def setup_model_directory(self):
        """设置模型保存目录"""
        if self.model_dir is None:
            # 默认保存到当前目录下的model文件夹
            self.model_dir = os.path.join(os.getcwd(), "model")

        # 创建目录（如果不存在）
        os.makedirs(self.model_dir, exist_ok=True)
        print(f"模型将保存到: {self.model_dir}")
        return self.model_dir

    def standardize_reanno_annotations(self, adata):
        """
        标准化reanno注释，将Neural crest；Neural.crest都定义为：Neural.crest
        """
        print("\n=== 标准化reanno注释 ===")

        if 'reanno' not in adata.obs.columns:
            print("警告: 未找到reanno列，跳过标准化")
            return adata

        # 获取原始的reanno值分布
        original_counts = adata.obs['reanno'].value_counts()
        print("原始reanno分布:")
        for cell_type, count in original_counts.items():
            print(f"  {cell_type}: {count}")

        # 将reanno列转换为字符串类型，避免Categorical问题
        reanno_series = adata.obs['reanno'].astype(str)
        print("将所有空格替换为点号...")
        reanno_series = reanno_series.str.replace(' ', '.', regex=False)

        # 标准化Neural crest相关注释
        reanno_series = reanno_series.replace({
            'Neural crest': 'Neural.crest',
            'Oligo': 'Oligodendrocyte'
        })

        reanno_series = reanno_series.str.replace(r'\.+', '.', regex=True)  # 将多个点号替换为单个点号
        reanno_series = reanno_series.str.strip('.')  # 去除首尾的点号
        # 更新reanno列
        adata.obs['reanno'] = reanno_series

        # 显示标准化后的分布
        standardized_counts = adata.obs['reanno'].value_counts()
        print("标准化后reanno分布:")
        for cell_type, count in standardized_counts.items():
            print(f"  {cell_type}: {count}")

        return adata

    def check_raw_counts_data(self, adata, dataset_name=""):
        """
        检查数据是否是原始counts数据

        参数:
        - adata: AnnData对象
        - dataset_name: 数据集名称用于标识

        返回:
        - is_raw_counts: 是否是原始counts数据的布尔值
        """
        print(f"\n=== 检查 {dataset_name} 的原始数据状态 ===")

        # 检查是否有layers['counts']
        has_counts_layer = 'counts' in adata.layers
        print(f"是否有counts层: {has_counts_layer}")

        if has_counts_layer:
            counts_data = adata.layers['counts']
            # 检查counts层的数据类型和值范围
            if hasattr(counts_data, 'toarray'):
                counts_data = counts_data.toarray()

            # 检查是否为整数
            is_integer = np.all(counts_data.astype(int) == counts_data)
            min_val = np.min(counts_data)
            max_val = np.max(counts_data)
            mean_val = np.mean(counts_data)

            print(f"counts层 - 是否为整数: {is_integer}")
            print(f"counts层 - 最小值: {min_val}, 最大值: {max_val}, 平均值: {mean_val:.2f}")

            # 判断是否为原始counts数据的标准
            is_raw_counts = is_integer and min_val >= 0 and mean_val > 0.1
            print(f"counts层是否为原始counts数据: {is_raw_counts}")

            return is_raw_counts

        else:
            print("没有counts层，检查X矩阵...")
            # 检查X矩阵
            x_data = adata.X
            if hasattr(x_data, 'toarray'):
                x_data = x_data.toarray()

            # 检查是否为整数
            is_integer = np.all(x_data.astype(int) == x_data)
            min_val = np.min(x_data)
            max_val = np.max(x_data)
            mean_val = np.mean(x_data)

            print(f"X矩阵 - 是否为整数: {is_integer}")
            print(f"X矩阵 - 最小值: {min_val}, 最大值: {max_val}, 平均值: {mean_val:.2f}")

            # 判断是否为原始counts数据的标准
            is_raw_counts = is_integer and min_val >= 0 and mean_val > 0.1
            print(f"X矩阵是否为原始counts数据: {is_raw_counts}")

            return is_raw_counts

    def fix_duplicate_index(self, adata, dataset_name=""):
        """彻底修复重复的索引"""
        print(f"检查 {dataset_name} 的索引唯一性...")

        # 检查是否有重复索引
        if adata.obs.index.duplicated().any():
            duplicate_count = adata.obs.index.duplicated().sum()
            print(f"发现 {duplicate_count} 个重复索引，正在修复...")

            # 获取重复的索引值
            duplicate_indices = adata.obs.index[adata.obs.index.duplicated()].unique()
            print(f"重复的索引示例: {list(duplicate_indices)[:5]}")

            # 创建新的唯一索引
            new_index = []
            count_dict = {}

            for original_idx in adata.obs.index:
                if original_idx not in count_dict:
                    count_dict[original_idx] = 0
                    new_index.append(original_idx)
                else:
                    count_dict[original_idx] += 1
                    new_index.append(f"{original_idx}_dup{count_dict[original_idx]}")

            # 更新索引
            adata.obs.index = new_index
            print("重复索引修复完成")
        else:
            print(f"{dataset_name} 索引唯一性检查通过")

        return adata

    def debug_index_issues(self, adata, name=""):
        """调试索引问题"""
        print(f"\n=== 调试 {name} 的索引 ===")
        print(f"索引类型: {type(adata.obs.index)}")
        print(f"索引长度: {len(adata.obs.index)}")
        print(f"唯一索引数量: {len(adata.obs.index.unique())}")
        print(f"是否有重复: {adata.obs.index.duplicated().any()}")

        if adata.obs.index.duplicated().any():
            duplicates = adata.obs.index[adata.obs.index.duplicated()].unique()
            print(f"重复的索引值 (前10个): {list(duplicates)[:10]}")

        return adata

    def clean_data_types(self, adata):
        """彻底清理数据类型，确保数值列是数值类型"""
        print("清理数据类型...")

        # 定义数值列
        numeric_columns = ['nCount_RNA', 'nFeature_RNA', 'percent.mt']

        for col in numeric_columns:
            if col in adata.obs.columns:
                print(f"处理列: {col}")
                print(f"  原始类型: {adata.obs[col].dtype}")
                print(f"  前5个值: {adata.obs[col].head().tolist()}")

                # 检查是否有非数值内容
                if adata.obs[col].dtype == 'object':
                    # 尝试转换为数值
                    try:
                        adata.obs[col] = pd.to_numeric(adata.obs[col], errors='coerce')
                        print(f"  转换为数值类型: {adata.obs[col].dtype}")
                    except Exception as e:
                        print(f"  转换失败: {e}")
                        # 如果转换失败，删除该列
                        print(f"  删除列: {col}")
                        adata.obs.drop(columns=[col], inplace=True)
                else:
                    print(f"  已经是数值类型: {adata.obs[col].dtype}")

        # 清理其他对象类型列
        for col in adata.obs.columns:
            if adata.obs[col].dtype == 'object':
                # 检查是否可以转换为分类
                unique_ratio = adata.obs[col].nunique() / len(adata.obs[col])
                if unique_ratio < 0.5:  # 如果唯一值比例小于50%，转为分类
                    adata.obs[col] = adata.obs[col].astype('category')
                    print(f"  将 {col} 转换为分类类型")
                else:
                    # 确保字符串类型
                    adata.obs[col] = adata.obs[col].astype(str)
                    print(f"  将 {col} 转换为字符串类型")

        return adata

    def load_and_validate_data(self, data_path):
        """加载和验证数据"""
        print("加载数据...")

        # 加载主要数据
        print(f"加载主要数据: {data_path}")
        adata = sc.read_h5ad(data_path)
        print(f"主要数据维度: {adata.n_obs} 细胞, {adata.n_vars} 基因")

        # 检查原始数据状态
        self.check_raw_counts_data(adata, "主要数据")

        # 调试主要数据索引
        adata = self.debug_index_issues(adata, "主要数据")
        adata = self.fix_duplicate_index(adata, "主要数据")
        adata = self.clean_data_types(adata)

        # 标准化主要数据的reanno注释
        adata = self.standardize_reanno_annotations(adata)

        return adata

    def preprocess_data(self, adata):
        """数据预处理"""
        print("数据预处理...")

        # 确保模型目录已设置
        self.setup_model_directory()

        # 使用临时副本进行高变基因选择，避免数值计算问题
        print("使用临时副本进行高变基因选择...")
        adata_tmp = adata.copy()

        # 对临时数据进行标准化和log转换
        sc.pp.normalize_total(adata_tmp, target_sum=1e4)
        sc.pp.log1p(adata_tmp)

        # 选择高变基因
        sc.pp.highly_variable_genes(
            adata_tmp,
            n_top_genes=self.n_top_genes,
            flavor="cell_ranger",
            batch_key="orig.ident",
            subset=False
        )

        # 使用高变基因筛选原始数据
        print(f"使用 {self.n_top_genes} 个高变基因筛选数据...")
        hvg_mask = adata_tmp.var['highly_variable']
        adata_hvg = adata[:, hvg_mask].copy()

        # 将高变基因信息添加到原始数据
        adata.var['highly_variable'] = hvg_mask

        # 保存预处理数据前再次清理数据类型
        adata = self.clean_data_types(adata)

        # 确保reanno注释已标准化
        adata = self.standardize_reanno_annotations(adata)

        # 准备用于scPoli训练的数据 - 手动处理稀疏矩阵
        print("处理稀疏矩阵...")
        if hasattr(adata_hvg.X, 'toarray'):
            # 如果是稀疏矩阵，转换为稠密矩阵
            adata_hvg.X = adata_hvg.X.toarray()
            print("已将稀疏矩阵转换为稠密矩阵")
        
        # 确保counts层存在
        if "counts" not in adata_hvg.layers:
            adata_hvg.layers["counts"] = adata_hvg.X.copy()
            print("创建了counts层")
        else:
            # 如果counts层存在但也是稀疏矩阵，同样需要转换
            if hasattr(adata_hvg.layers["counts"], 'toarray'):
                adata_hvg.layers["counts"] = adata_hvg.layers["counts"].toarray()
            adata_hvg.X = adata_hvg.layers["counts"].copy()

        print(f"高变基因数据维度: {adata_hvg.shape}")
        print(f"高变基因数量: {adata_hvg.n_vars}")

        return adata_hvg

    def train_scpoli_model(self, adata):
        """训练scPoli模型"""
        # 设置模型目录
        model_dir = self.setup_model_directory()
        model_subdir = os.path.join(model_dir, f'scpoli_model_{self.cell_type_key}_hvg{self.n_top_genes}')

        print("\n=== scPoli 模型训练 ===")

        # 设置参数
        condition_key = "orig.ident"
        cell_type_key = self.cell_type_key

        print(f"条件变量: {condition_key}")
        print(f"细胞类型变量: {cell_type_key}")
        print(f"高变基因数量: {self.n_top_genes}")
        print(f"条件分布:\n{adata.obs[condition_key].value_counts()}")
        print(f"细胞类型分布:\n{adata.obs[cell_type_key].value_counts()}")

        early_stopping_kwargs = {
            "early_stopping_metric": "val_prototype_loss",
            "mode": "min",
            "threshold": 0,
            "patience": 20,
            "reduce_lr": True,
            "lr_patience": 13,
            "lr_factor": 0.1,
        }

        # 训练scPoli模型
        print("训练scPoli模型...")

        scpoli_model = scPoli(
            adata=adata,
            condition_keys=condition_key,
            cell_type_keys=cell_type_key,
            embedding_dims=128,
            recon_loss='nb'
        )

        scpoli_model.train(
            n_epochs=200,
            pretraining_epochs=40,
            early_stopping_kwargs=early_stopping_kwargs,
            eta=5
        )

        return scpoli_model, model_subdir

    def save_model_and_results(self, scpoli_model, model_dir, adata):
        """保存模型和结果"""
        # 修复数据类型
        adata = self.clean_data_types(adata)

        # 确保reanno注释已标准化
        adata = self.standardize_reanno_annotations(adata)

        # 保存模型
        try:
            print("保存模型...")
            scpoli_model.save(model_dir, overwrite=True, save_anndata=True)
            print(f"✅ 模型已保存至 {model_dir}")
        except Exception as e:
            print(f"❌ 保存错误: {e}")

        return scpoli_model

    def get_latent_representation(self, scpoli_model, adata):
        """获取潜在表示"""
        print("获取scPoli潜在表示...")
        scpoli_model.model.eval()

        if hasattr(adata.X, 'toarray'):
            adata.X = adata.X.toarray()

        data_latent_source = scpoli_model.get_latent(adata, mean=True)
        adata.obsm["scPoli"] = data_latent_source

        return adata

    def visualize_results(self, adata):
        """可视化结果"""
        print("基于scPoli的聚类和可视化...")

        # 确保模型目录已设置
        self.setup_model_directory()

        # 计算邻居和UMAP
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.neighbors(adata, use_rep="scPoli", random_state=self.seed)
        sc.tl.umap(adata, random_state=self.seed)

        # 绘制UMAP图
        plt.figure(figsize=(6, 5))
        sc.pl.umap(adata, color="orig.ident", frameon=False, show=False,
                title=f"orig.ident - {self.n_top_genes} HVGs")
        plt.tight_layout()

        # 保存到模型目录
        plot_path = os.path.join(self.model_dir, f"scpoli_orig_ident_{self.cell_type_key}_hvg{self.n_top_genes}.pdf")
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        plt.close()

        plt.figure(figsize=(6, 5))
        sc.pl.umap(adata, color=self.cell_type_key, frameon=False, show=False,
                title=f"{self.cell_type_key} - {self.n_top_genes} HVGs")
        plt.tight_layout()

        plot_path = os.path.join(self.model_dir, f"scpoli_{self.cell_type_key}_hvg{self.n_top_genes}.pdf")
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        plt.close()

        # 保存整合后的数据到模型目录
        output_path = os.path.join(self.model_dir, f"scpoli_integrated_{self.cell_type_key}_hvg{self.n_top_genes}.h5ad")

        # 确保数据类型正确
        adata = self.clean_data_types(adata)
        # 确保reanno注释已标准化
        adata = self.standardize_reanno_annotations(adata)

        try:
            adata.write_h5ad(output_path)
            print(f"整合数据已保存至: {output_path}")
        except Exception as e:
            print(f"保存整合数据时出错: {e}")
            # 尝试保存到当前目录
            fallback_path = f'scpoli_integrated_{self.cell_type_key}_hvg{self.n_top_genes}.h5ad'
            try:
                adata.write_h5ad(fallback_path)
                print(f"整合数据已保存到当前目录: {fallback_path}")
            except Exception as e2:
                print(f"整合数据保存失败: {e2}")

    def run_full_analysis(self):
        """运行完整分析流程"""
        # 数据路径
        data_path = "/storage2/liuxiaodongLab/jiangjing/Projects/XueyingFan/PD_XueyingFan/20251106_fetal_adult_merge_name_scdevelopment/output/filter_defined_number_combined_sampled.h5ad"

        try:
            # 1. 加载数据
            adata = self.load_and_validate_data(data_path)

            # 2. 数据预处理
            adata_hvg = self.preprocess_data(adata)

            # 3. 训练scPoli模型
            scpoli_model, model_dir = self.train_scpoli_model(adata_hvg)

            # 4. 保存模型
            scpoli_model = self.save_model_and_results(scpoli_model, model_dir, adata_hvg)

            # 5. 获取潜在表示
            adata_hvg = self.get_latent_representation(scpoli_model, adata_hvg)

            # 6. 可视化结果
            self.visualize_results(adata_hvg)

            print("✅ 分析完成！")

        except Exception as e:
            print(f"❌ 分析过程中出现错误: {e}")
            raise

def main():
    """主函数，支持命令行参数"""
    parser = argparse.ArgumentParser(description='scPoli整合分析')
    parser.add_argument('--cell_type_key', type=str, choices=['reanno', 'lineage'],
                       default='reanno', help='细胞类型键: reanno 或 lineage (默认: reanno)')
    parser.add_argument('--n_top_genes', type=int, choices=[2000, 4000],
                       default=2000, help='高变基因数量: 2000 或 4000 (默认: 2000)')
    parser.add_argument('--seed', type=int, default=42, help='随机种子 (默认: 42)')
    parser.add_argument('--model_dir', type=str, default=None,
                       help='模型保存目录 (默认: 当前目录下的model文件夹)')

    args = parser.parse_args()

    print(f"开始scPoli分析")
    print(f"细胞类型键: {args.cell_type_key}")
    print(f"高变基因数: {args.n_top_genes}")
    print(f"随机种子: {args.seed}")
    print(f"模型目录: {args.model_dir if args.model_dir else '当前目录下的model文件夹'}")

    analyzer = scPoliIntegration(
        seed=args.seed,
        cell_type_key=args.cell_type_key,
        n_top_genes=args.n_top_genes,
        model_dir=args.model_dir
    )
    analyzer.run_full_analysis()

# 运行分析
if __name__ == "__main__":
    main()
