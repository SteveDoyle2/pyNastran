class KeyMap:
    def __init__(self, subtitle: str,
                 label: str,
                 superelement_adaptivity_index: str,
                 pval_step: str):
        self.subtitle = subtitle
        self.label = label # .replace(' '*20, ' ')
        self.superelement_adaptivity_index = superelement_adaptivity_index
        self.pval_step = pval_step

KeysMap = dict[str, KeyMap]
