$(eval $(begin-module))

################################################################
# unit definitions
################################################################

module_units_h += utility constants tprme
module_units_cpp-h := relativecm_rme
# module_units_f :=

module_programs_cpp += relativecm-gen
module_programs_cpp_test := relativecm_rme_test
# module_programs_f :=
# module_generated :=

################################################################
# library creation flag
################################################################

$(eval $(library))

$(eval $(end-module))
