import pandas as pd  # type: ignore
#import xlwings as xl

def debug_output(objects, to_excel=True):
    variables=[attr for attr in dir(objects[0]) if not callable(attr) and not attr.startswith("__")]

    df = pd.DataFrame([[getattr(i,j) for j in variables] for i in objects], columns = variables)
    if to_excel:
        df.to_clipboard()
        #to_excel('xdb_objects.xlsx')
        #sht = xl.Book().sheets[0]
        #sht.range('A1').value=df
    return df

