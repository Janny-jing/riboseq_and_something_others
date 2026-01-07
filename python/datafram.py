import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
import numpy as np

def should_remove(text, remove_strings):
    """检查文本是否包含需要移除的字符串"""
    text_lower = str(text).lower()
    return any(remove_string.lower() in text_lower for remove_string in remove_strings)

def clean_filename(filename, remove_strings):
    """去掉文件名中的指定字符串，并去除扩展名"""
    base_name = os.path.splitext(filename)[0]  # 去掉扩展名
    if should_remove(base_name, remove_strings):
        return None
    return base_name

def clean_column_name(col_name):
    """清理列名，去除'_new_sum'部分"""
    return col_name.replace('_new_sum', '').strip()

def combine_xlsx_files(input_folder, output_file, heatmap_file, remove_strings, font_size, color_map):
    # 创建一个空的 DataFrame 来存储结果
    result_df = pd.DataFrame()

    # 遍历文件夹中的所有文件
    for filename in os.listdir(input_folder):
        if filename.endswith('.xlsx'):
            cleaned_filename = clean_filename(filename, remove_strings)
            if cleaned_filename is None:
                continue
            
            # 构建完整的文件路径
            file_path = os.path.join(input_folder, filename)
            
            # 读取 Excel 文件
            df = pd.read_excel(file_path, engine='openpyxl')
            
            # 提取第 5 列和第 6 列数据
            column_5 = df.iloc[:, 4]  # 第 5 列索引为 4
            column_6 = df.iloc[:, 5]  # 第 6 列索引为 5
            
            # 筛选行名不包含需要移除的字符串的行
            filtered_df = df[~column_5.apply(lambda x: should_remove(x, remove_strings))]
            column_5_filtered = filtered_df.iloc[:, 4]
            column_6_filtered = filtered_df.iloc[:, 5]

            # 将数据合并到一个临时 DataFrame 中，并使用清理后的文件名作为列名
            temp_df = pd.DataFrame({
                'Key': column_5_filtered,
                cleaned_filename: column_6_filtered
            })
            
            # 如果是第一个文件，则直接赋值给 result_df
            if result_df.empty:
                result_df = temp_df
            else:
                # 合并当前文件的数据到 result_df 中，基于 'Key' 列
                result_df = pd.merge(result_df, temp_df, on='Key', how='outer')

    # 确保 'Key' 列不被转换为数字索引
    result_df.set_index('Key', inplace=True)

    # 过滤掉包含需要移除字符串的列
    columns_to_keep = [col for col in result_df.columns if not should_remove(col, remove_strings)]
    result_df = result_df[columns_to_keep]

    # 清理列名中的 '_new_sum'
    result_df.columns = [clean_column_name(col) for col in result_df.columns]

    # 处理缺失值：将 NaN 替换为 0 或其他合适的值
    result_df.fillna(0, inplace=True)

    # 保存结果到新的 Excel 文件
    result_df.to_excel(output_file)
    print(f"Data has been saved to {output_file}")

    # 设置全局字体大小
    plt.rcParams.update({'font.size': font_size})

    # 绘制热图并保存
    plt.figure(figsize=(50,40))
    
    # 检查并转换数据类型为浮点数
    heatmap_data = result_df.astype(float)

    # 使用 mask 来隐藏缺失值或特定条件下的值
    mask = heatmap_data.isnull()
    
    # 绘制热图
    sns.heatmap(heatmap_data, cmap=color_map, annot=False, mask=mask, cbar_kws={'label': 'Value'})
    
    # 设置标题和标签字体大小
    plt.title('Heatmap of Combined Data', fontsize=font_size + 2)
    plt.xticks(fontsize=font_size)
    plt.yticks(fontsize=font_size)
    
    plt.tight_layout()
    
    # 保存热图到指定文件
    plt.savefig(heatmap_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Heatmap saved to {heatmap_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Combine multiple xlsx files into one for heatmap.")
    parser.add_argument("input_folder", help="Path to the folder containing the xlsx files.")
    parser.add_argument("output_file", help="Name of the output Excel file (with .xlsx extension).")
    parser.add_argument("heatmap_file", help="Name of the output heatmap image file (with .png extension).")
    parser.add_argument("--remove", nargs='+', default=[], 
                        help="Strings to remove from row and column names.")
    parser.add_argument("--font-size", type=int, default=12, 
                        help="Font size for labels in the heatmap.")
    parser.add_argument("--color-map", type=str, default='viridis', 
                        help="Color map for the heatmap (e.g., viridis, plasma, coolwarm).")
    
    args = parser.parse_args()

    combine_xlsx_files(args.input_folder, args.output_file, args.heatmap_file, args.remove, args.font_size, args.color_map)
