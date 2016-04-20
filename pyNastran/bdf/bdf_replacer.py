from __future__ import (nested_scopes, absolute_import,
                        print_function, unicode_literals)
from six import PY2
from pyNastran.bdf.field_writer_8 import print_card_8
from pyNastran.bdf.bdf_interface.bdf_card import BDFCard
from pyNastran.bdf.cards.utils import wipe_empty_fields
from pyNastran.bdf.bdf_interface.assign_type import interpret_value
from pyNastran.bdf.bdf import BDF, to_fields


class BDFReplacer(BDF):
    """
    This class demonstates OpenMDAO find-replace streaming coupled with some
    Python black magic (see _read_bulk_data_deck).
    """
    def __init__(self, bdf_out_filename, debug=True, log=None):
        BDF.__init__(self, debug=debug, log=log)
        self.bdf_out_filename = bdf_out_filename

    def is_reject(self, card_name):
        return False

    def _start_writing(self):
        if PY2:
            self.bdf_out_file = open(self.bdf_out_filename, 'wb')
        else:
            self.bdf_out_file = open(self.bdf_out_filename, 'w')
        self._write_header(self.bdf_out_file)

    def _read_bulk_data_deck(self):
        """
        Hacks the _read_bulk_data_deck method to do some pre/post processing,
        but still calls the orignal function.

        ..see::  BDF._read_bulk_data_deck()
        """
        self._start_writing()
        BDF._read_bulk_data_deck(self)
        self._finish_writing()

    def _finish_writing(self):
        self.bdf_out_file.write('ENDDATA\n')
        self.bdf_out_file.close()

    def add_card(self, card_lines, card_name, comment='', is_list=True):
        """
        Adds a card object to the BDF object using a streaming approach.

        Parameters
        ----------
        card_lines : List[str]
            the list of the card fields

         >>> ['GRID,1,2',]  # (is_list = False)
         >>> ['GRID',1,2,]  # (is_list = True; default)

        card_name : str
            the card_name -> 'GRID'
        comment : str; default=True
            an optional the comment for the card
        is_list : bool; default=True
            changes card_lines from a list of lines to
            a list of fields

        Returns
        -------
        card_object : BDFCard()
            the card object representation of card

        .. note:: this is a very useful method for interfacing with the code
        .. note:: the cardObject is not a card-type object...so not a GRID
                  card or CQUAD4 object.  It's a BDFCard Object.  However,
                  you know the type (assuming a GRID), so just call the
                  *mesh.Node(nid)* to get the Node object that was just
                  created.
        .. warning:: cardObject is not returned
        """
        if comment:
            self.bdf_out_file.write(comment)

        if card_name in ['DEQATN']:
            card_obj = card_lines
            card = card_lines
            #print("card =", '\n'.join(card_lines))
            self.bdf_out_file.write('\n'.join(card_lines))
        else:
            if is_list:
                fields = card_lines
            else:
                fields = to_fields(card_lines, card_name)

            # apply OPENMDAO syntax
            if self._is_dynamic_syntax:
                fields = [self._parse_dynamic_syntax(field) if '%' in
                          field[0:1] else field for field in fields]

            card = wipe_empty_fields([interpret_value(field, fields) if field is not None
                                      else None for field in fields])
            #print(card)
            card_obj = BDFCard(card)
            self.bdf_out_file.write(print_card_8(card_obj))

        #for reject_card in self.reject_cards:
            #try:
            #print('comment =', comment)
            #except RuntimeError:
                #for field in reject_card:
                    #if field is not None and '=' in field:
                        #raise SyntaxError('cannot reject equal signed '
                        #              'cards\ncard=%s\n' % reject_card)
                #raise

        #self.reject_cards.append(card)
        #print('rejecting processed auto=rejected %s' % card)
        return card_obj
