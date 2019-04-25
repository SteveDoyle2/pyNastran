from pyNastran.op2.op2 import OP2


model = OP2()
methods_dict = {
  b'OEFITSTN' : [model._table_passer, model._table_passer],
  #b'OEFITSTN' : False,  # another option
}
model.set_additional_result_tables_to_read(methods_dict)
model.read_op2('../../../models/elements/static_elements.op2')
model.export_to_hdf5('static_elements.h5')
