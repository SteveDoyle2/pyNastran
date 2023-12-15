class KeyMap:
    def __init__(self, subtitle: str,
                 label: str,
                 superelement_adaptivity_index: str,
                 pval_step: str):
        self.subtitle = subtitle

        end_label = label[-20:]  # TODO: is 20 the right number?
        if end_label.strip():
            slabel = end_label.split()
            if len(slabel) == 2 and slabel[0] == 'SUBCASE' and slabel[1].isdigit():
                label = label[:-20].strip()
                # in:  'TIP CENTER LOAD                                                                                        SUBCASE 1'
                # out: 'TIP CENTER LOAD'
        self.label = label # .replace(' '*20, ' ')
        self.superelement_adaptivity_index = superelement_adaptivity_index
        self.pval_step = pval_step

    def __repr__(self) -> str:
        msg = 'KeyMap:\n'
        if self.subtitle:
            msg += f'  subtitle: {subtitle.label!r}\n'
        if self.label:
            msg += f'  label: {self.label!r}\n'
        if self.superelement_adaptivity_index:
            msg += f'  superelement_adaptivity_index: {self.superelement_adaptivity_index!r}\n'
        if self.pval_step:
            msg += f'  pval_step: {self.pval_step!r}'

        if msg == 'KeyMap:\n':
            msg = 'KeyMap(...)'
        return msg

KeysMap = dict[str, KeyMap]

# ints: subcase, analysis_code, sort_code, count, ogs
# strs: superelement_adaptivity_index, pval_step
NastranKey = tuple[int, int, int, int, int,
                   str, str]
