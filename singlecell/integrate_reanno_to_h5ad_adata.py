import pandas as pd
import os
import glob
import argparse
import sys
import scanpy as sc

def integrate_reanno_predictions(csv_folder, original_adata_path, output_adata_path):
    """
    整合文件夹中所有CSV文件的reanno_pred列到原单细胞数据的obs中
    
    参数:
    csv_folder: 包含CSV文件的文件夹路径
    original_adata_path: 原单细胞数据文件路径(.h5ad)
    output_adata_path: 输出单细胞数据文件路径(.h5ad)
    """
    
    # 检查输入文件是否存在
    if not os.path.exists(original_adata_path):
        print(f"错误: 原单细胞数据文件不存在: {original_adata_path}")
        sys.exit(1)
    
    if not os.path.exists(csv_folder):
        print(f"错误: CSV文件夹不存在: {csv_folder}")
        sys.exit(1)
    
    # 获取所有CSV文件
    csv_files = glob.glob(os.path.join(csv_folder, "*.csv"))
    if len(csv_files) == 0:
        print(f"错误: 在文件夹 {csv_folder} 中未找到任何CSV文件")
        sys.exit(1)
    
    print(f"找到 {len(csv_files)} 个CSV文件")
    
    # 读取原单细胞数据
    print(f"读取原单细胞数据: {original_adata_path}")
    adata = sc.read(original_adata_path)
    print(f"原数据包含 {adata.n_obs} 个细胞, {adata.n_vars} 个基因")
    
    # 用于存储所有reanno_pred数据
    all_reanno_data = {}
    
    # 处理每个CSV文件
    for i, csv_file in enumerate(csv_files):
        filename = os.path.basename(csv_file)
        print(f"处理文件 {i+1}/{len(csv_files)}: {filename}")
        
        try:
            # 读取CSV文件
            df = pd.read_csv(csv_file, index_col=0)
            
            # 检查索引是否匹配
            common_cells = df.index.intersection(adata.obs.index)
            if len(common_cells) == 0:
                print(f"  警告: {filename} 的细胞索引与原数据完全不匹配")
                continue
            elif len(common_cells) < len(adata.obs.index):
                print(f"  注意: {filename} 只包含 {len(common_cells)}/{len(adata.obs.index)} 个共同细胞")
            
            # 提取reanno_pred列
            if 'reanno_pred' in df.columns:
                # 使用文件名作为新列名（去除扩展名）
                col_name = f"reanno_pred_{os.path.splitext(filename)[0]}"
                
                # 重新索引以匹配原数据
                reanno_series = df['reanno_pred'].reindex(adata.obs.index)
                
                # 添加到原数据的obs中
                adata.obs[col_name] = reanno_series
                print(f"  成功添加列: {col_name}")
                
                # 统计非空值
                non_null_count = reanno_series.notna().sum()
                print(f"  有效预测值: {non_null_count}/{len(adata.obs)} 个细胞")
                
            else:
                print(f"  警告: {filename} 中未找到 reanno_pred 列")
                
        except Exception as e:
            print(f"  错误: 处理 {filename} 时出错 - {e}")
            continue
    
    # 保存新的h5ad文件
    print(f"\n保存整合后的单细胞数据到: {output_adata_path}")
    adata.write(output_adata_path)
    
    # 显示统计信息
    print(f"\n整合完成！")
    print(f"原单细胞数据已添加 {len([col for col in adata.obs.columns if col.startswith('reanno_pred_')])} 个reanno_pred列")
    print(f"新列名: {[col for col in adata.obs.columns if col.startswith('reanno_pred_')]}")
    
    # 显示新添加的列的前几行
    reanno_cols = [col for col in adata.obs.columns if col.startswith('reanno_pred_')]
    if reanno_cols:
        print(f"\n新添加的列数据预览:")
        print(adata.obs[reanno_cols].head())
    
    return adata

def main():
    parser = argparse.ArgumentParser(description='整合多个CSV文件中的reanno_pred列到单细胞数据中')
    parser.add_argument('--csv_folder', '-c', required=True, help='包含CSV文件的文件夹路径')
    parser.add_argument('--original_adata', '-o', required=True, help='原单细胞数据文件路径(.h5ad)')
    parser.add_argument('--output_adata', '-out', required=True, help='输出单细胞数据文件路径(.h5ad)')
    
    args = parser.parse_args()
    
    print("开始整合reanno_pred数据...")
    print(f"CSV文件夹: {args.csv_folder}")
    print(f"原单细胞数据: {args.original_adata}")
    print(f"输出文件: {args.output_adata}")
    print("-" * 50)
    
    integrate_reanno_predictions(args.csv_folder, args.original_adata, args.output_adata)

if __name__ == "__main__":
    main()
