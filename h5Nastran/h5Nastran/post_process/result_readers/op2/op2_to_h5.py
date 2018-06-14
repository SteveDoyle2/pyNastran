from pyNastran.op2.op2 import OP2


model = OP2()
methods_dict = {
  b'OEFITSTN' : [model._table_passer, model._table_passer],
  #b'OEFITSTN' : False,  # another option
}
model.set_additional_result_tables_to_read(methods_dict)
# model.read_op2(r'\\na\shares\Structure\GAR1-STR-SA\05_FEA\HJ1_IFEM03\Output_C7\10144_Rev_E\Output\OP2\HA420_IFEM03_C7UL_BAM3_01.op2')
# model.export_to_hdf5('HA420_IFEM03_C7UL_BAM3_01.h5')
model.read_op2('../../../models/elements/static_elements.op2')
model.export_to_hdf5('static_elements.h5')
