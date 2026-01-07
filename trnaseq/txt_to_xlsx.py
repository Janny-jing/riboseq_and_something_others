import argparse
import pandas as pd

def txt_to_xlsx(input_file, output_file):
    # 读取文本文件
    df = pd.read_csv(input_file, sep='\t')
    
    # 将 DataFrame 写入 Excel 文件
    df.to_excel(output_file, index=False)

if __name__ == "__main__":
    # 创建 ArgumentParser 对象
    parser = argparse.ArgumentParser(description="Convert a .txt file to an .xlsx file.")
    
    # 添加命令行参数
    parser.add_argument("input", help="The input .txt file")
    parser.add_argument("output", help="The output .xlsx file")
    
    # 解析命令行参数
    args = parser.parse_args()
    
    # 调用函数
    txt_to_xlsx(args.input, args.output)
