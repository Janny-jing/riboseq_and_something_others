# 指定要处理的文件夹路径和输出文件夹路径
$inputFolderPath = "F:\work\trnaid_1000\"
$outputFolderPath = "F:\work\trnaid_1000_out\"

# 获取文件夹中所有的文件
$fileList = Get-ChildItem $inputFolderPath -Filter "*.txt"

# 遍历每个文件
foreach ($file in $fileList) {
    # 读取文件内容
    $fileContent = Get-Content $file.FullName
    
    # 处理每一行
    $processedLines = $fileContent | ForEach-Object {
        # 替换分隔符 ":"
        $line = $_.Trim() -replace ':', "`n" -replace ',', "`n" -replace ' ', ''

        # 将替换后的字符串按换行分割成数组
        $parts = $line.Split("`n")
        
        # 返回处理后的数组
        $parts
    }
    
    # 删除重复行
    $uniqueLines = $processedLines | Select-Object -Unique
    
    # 构建输出文件路径
    $outputFilePath = Join-Path $outputFolderPath $file.Name
    
    # 输出结果到新的文件
    $uniqueLines | Out-File $outputFilePath -Encoding UTF8
}