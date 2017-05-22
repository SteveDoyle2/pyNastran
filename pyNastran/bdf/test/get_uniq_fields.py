"""
defines BDFUniqueCard
"""
from __future__ import print_function
import sys
import logging

from pyNastran.bdf.bdf_interface.assign_type import interpret_value
from pyNastran.bdf.bdf import BDF, to_fields, BDFCard
#from pyNastran.bdf.cards.utils import wipe_empty_fields


class BDFUniqueCard(BDF):
    """
    class that creates a list of files that use a different form of the
    cards to try to focus in on problems better.  For example, a GRID
    can be of the form:

    GRID,1
    GRID,2,, 0.1, 0.2, 0.3
    GRID,3,1,0.1, 0.2, 0.3

    This code will find the first instance of each
    """
    def __init__(self, log, card_fingerprint_set):
        """
        Initializes the BDFUniqueCard object

        Parameters
        ----------
        #debug : bool/None; default=True
            #used to set the logger if no logger is passed in
                #True:  logs debug/info/error messages
                #False: logs info/error messages
                #None:  logs error messages
        log : logging module object / None
            if log is set, debug is ignored and uses the
            settings the logging object has
        card_fingerprint_set : set([str, str, ...])
            a set to store data
        """
        self.card_set = card_fingerprint_set
        #self.new_bdf = new_bdf
        BDF.__init__(self, log=log)

    #def _parse_executive_control_deck(self):
      #pass

    def cross_reference(self, *args, **kwargs):
        """fakes cross_reference"""
        pass

    def add_card(self, card_lines, card_name, comment=''):
        card = to_fields(card_lines, card_name)
        card = [interpret_value(val) for val in card]
        card = BDFCard(card)

        #if card_name == "LOAD":
        rec = []
        for item in card:
            if isinstance(item, int):
                rec.append(1)
            elif isinstance(item, float):
                rec.append(1.0)
            else:
                rec.append(item)
            #print(type(c), c)
        #print("*", rec)
        rec_str = str(rec)
        if not rec_str in self.card_set:
            #self.new_bdf.add_card(card, card_name, icard, old_card_obj)
            #print(card)
            self.card_set.add(rec_str)

    def is_reject(self, card_name):
        """Can the card be read"""
        if card_name.startswith('='):
            return False
        elif not card_name in self.cards_to_read:
            if card_name:
                if card_name not in self.reject_count:
                    self.log.info("reject_card_name = %r" % card_name)
                    self.reject_count[card_name] = 0
                self.reject_count[card_name] += 1
        return False

    def write_bdf(self, *args, **kwargs):
        """fakes write_bdf"""
        with open('cards.out.bdf', 'w') as bdf_file:
            for card in self.card_set:
                bdf_file.write(str(card) + '\n')
        print("please see cards.out.bdf")

def get_unique_bdf_cards(infilenames, logger):
    """runs BDFUniqueCard"""
    print("enter list of filenames at the command line...")
    card_fingerprint_set = set()
    assert len(infilenames) > 0, infilenames
    for infilename in infilenames:
        print("running infilename=%r" % infilename)
        try:
            infilename = infilename.strip()
            sys.stdout.flush()
            model = BDFUniqueCard(
                log=logger,
                card_fingerprint_set=card_fingerprint_set)
            model.read_bdf(infilename)
        except Exception as e:
            logger.error(e)
    model.write_bdf()
    print("done...")

def main():
    """runs BDFUniqueCard"""
    infilenames = sys.argv[1:]

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

    get_unique_bdf_cards(infilenames, logger)

if __name__ == '__main__':  # pragma: no cover
    main()
