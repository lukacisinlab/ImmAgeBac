import pdfplumber
import pandas as pd

pdf_file = "paper_exosomes.pdf"
tables=[]

table_settings = {
    "vertical_strategy": "explicit",
    "horizontal_strategy": "lines",
    "text_y_tolerance": 7,
    "explicit_vertical_lines": [34, 110, 278, 380, 612 - 34] #1 - [34, 120, 304, 394, 612 - 34] #2 - [34, 110, 278, 381, 612 - 34]
}



with pdfplumber.open(pdf_file) as pdf:
    pages = pdf.pages
    table = []
    for i,pg in enumerate(pages):
        if 32 >= i >= 27:  #26 20 #10 19
            page_table = pages[i].extract_table(table_settings=table_settings)
            table += page_table
    df = pd.DataFrame(table)
    with pd.option_context('display.max_rows', None,
                       'display.max_columns', None,
                       'display.precision', 3,
                       ):
        print(df)
    df.to_csv("table_test.csv")