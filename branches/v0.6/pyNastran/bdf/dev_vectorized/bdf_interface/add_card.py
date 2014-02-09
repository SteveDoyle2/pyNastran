class AddCard(object):
    def __init__(self):
        pass

    def add_PARAM(self, param, allowOverwrites=False):
        key = param.key
        if key in self.params and not allowOverwrites:
            if not param.isSameCard(self.params[key]):
                #assert param.key not in self.params,'key=%s param=%s oldPARAM=%s' %(key,param,self.params[key])
                self.log.warning('key=%s param=%s oldPARAM=%s' %
                                (key, param, self.params[key]))
                self.params[key] = param
        else:
            self.params[key] = param
