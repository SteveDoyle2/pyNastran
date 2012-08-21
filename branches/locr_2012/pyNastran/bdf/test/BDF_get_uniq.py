#!/usr/bin/python
import sys
import pyNastran.bdf.bdf
import logging

logger = logging.getLogger("bdfuniq")
logger.setLevel(logging.DEBUG)

log_formater = logging.Formatter(
    fmt="%(name)-10s: %(levelname)-5s: %(message)s")

log_handler_stderr = logging.StreamHandler()
log_handler_stderr.setFormatter(log_formater)
logger.addHandler(log_handler_stderr)

log_handler_file = logging.FileHandler("bdf_get_uniq.log", "a")
log_handler_file.setFormatter(log_formater)
logger.addHandler(log_handler_file)

includeDir = None


class BDFuniqCard(pyNastran.bdf.bdf.BDF):
    def __init__(self, log, fingset):
        self.card_set = fingset
        #self.newBDF = newBDF
        pyNastran.bdf.bdf.BDF.__init__(self, log=log)

    #def _parse_executive_control_deck(self):
    #   pass

    def cross_reference(self, xref):
        pass

    def add_card(self, card, cardName, iCard=0, old_card_obj=None):
        #if cardName == "LOAD":
        rec = []
        for item in card:
            item_type = type(item)
            if item_type == int:
                item = 1
                #print "int"
            elif item_type == float:
                item = 1.0
                #print "float"
            rec.append(item)
            #print type(c), c
        #print "*", rec
        rec_str = str(rec)
        if not rec_str in self.card_set:
            #self.newBDF.add_card(card, cardName, iCard, old_card_obj)
            print card
            self.card_set.add(rec_str)

    def _is_reject(self, cardName):
        """Can the card be read"""
        #cardName = self._get_card_name(card)
        if cardName.startswith('='):
            return False
        elif not cardName in self.cardsToRead:
            if cardName:
                if cardName not in self.rejectCount:
                    self.log.info("RejectCardName = |%s|" % (cardName))
                    self.rejectCount[cardName] = 0
                self.rejectCount[cardName] += 1
        return False


#newBDF = pyNastran.bdf.bdf.BDF()

card_fingerprint_set = set()

for infilename in sys.stdin.readlines():
    try:
        infilename = infilename.strip()
        BDFuniqCard(log=logger,
                    fingset=card_fingerprint_set).readBDF(infilename)
    except BaseException as e:
        logger.error(e)

#newBDF.writeBDF()

print "done..."