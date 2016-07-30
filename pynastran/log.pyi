def make_log(display=False):
    def stderr_logging(self, typ, msg:str):
    def __init__(self, level='debug':str, encoding='utf-8':str, log_func=None):
    def properties(self):
    def debug(self, msg:str):
    def msg_typ(self, typ, msg:str):
    def simple_msg(self, msg, typ=None:str):
    def info(self, msg:str):
    def warning(self, msg:str):
    def error(self, msg:str):
    def exception(self, msg:str):
    def critical(self, msg:str):
def get_logger(log=None, level='debug':str, encoding='utf-8':str):
def get_logger2(log=None, debug=True:bool, encoding='utf-8':str):
