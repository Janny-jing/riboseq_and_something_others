import pandas as pd
import argparse

def txt_to_xlsx(txt_file_path, xlsx_file_path):
    # 读取txt文件内容
    with open(txt_file_path, 'r') as file:
        lines = file.readlines()

    all_data = []
    group_data = {'Sequence': [], 'Count': []}
    
    for line in lines:
        stripped_line = line.strip()
        if stripped_line.startswith('Group'):
            if group_data['Sequence']:
                df = pd.DataFrame(group_data)
                all_data.append(df)
                group_data = {'Sequence': [], 'Count': []}
        else:
            parts = stripped_line.split()
            if len(parts) == 2:
                group_data['Sequence'].append(parts[0])
                group_data['Count'].append(int(parts[1]))

    # 添加最后一组数据
    if group_data['Sequence']:
        df = pd.DataFrame(group_data)
        all_data.append(df)

    # 组合所有数据框，每组占用两列，并命名列名以便区分不同组
    combined_df = pd.concat([df.add_prefix(f'Group{idx+1}_') for idx, df in enumerate(all_data)], axis=1)

    # 写入Excel文件
    combined_df.to_excel(xlsx_file_path, index=False)

if __name__ == "__main__":
    # 设置命令行参数解析
    parser = argparse.ArgumentParser(description='Convert a TXT file to an XLSX file.')
    parser.add_argument('input', help='Path to the input TXT file')
    parser.add_argument('output', help='Path to the output XLSX file')

    args = parser.parse_args()

    # 调用转换函数
    txt_to_xlsx(args.input, args.output)

    print(f"Conversion completed. Output saved to {args.output}")
