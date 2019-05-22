#import os
#import struct
from pyNastran.dev.xdb.enums import SolutionType



class XDB_obj:
    def __init__(self, table_name, ints, debug=False, log=None):
        self.raw_control = ints[:10]
        self.raw_full = ints
        self.raw_attributes = ints[10:]

        self.name = table_name
        self.project_no = ints[0]
        self.superelement_id = ints[1]
        self.path_no = ints[2]
        self.subcase_id = ints[3]
        self.set_id = ints[4]
        self.solution_code = SolutionType(ints[5])
        self.design_cycle = ints[6]
        self.iteration_cycle = ints[7]
        self.symmetry_segment = ints[8]
        # ints[9] - Intentionally left undefined

        #Second line
        self.internal_block_number = ints[10]
        self.data_block_format = ints[11]
        self.block_number_of_primary_map = ints[12] # Block number of the primary map for the object
        self.block_number_of_first_data_area = ints[13] # Block number of the first data area of the object
        self.no_of_words_per_entry = ints[14] # Number of words per entry, keyed objects only
        self.max_no_of_records_per_block = ints[15] # Maximum number of records per block, keyed objects only

        #Third line
        self.min_key_value = ints[16] # Minimum key value, keyed objects only
        self.entries_no = ints[17] # Either number of entries for keyed objects or number of words for sequential objects
        self.data_blocks_no = ints[18] # Number of data blocks
        self.secondary_index_blocks_no = ints[19] # Number of secondary index blocks

        self.max_key_value = ints[20] # Maximum key value, keyed objects only
        self.last_block_no = ints[21] # Block number of the last data block
