$(eval $(begin-module))

################################################################
# unit definitions
################################################################

module_units_h += chime constants tprme
module_units_cpp-h := relative_rme relativecm_rme
# module_units_f :=

module_programs_cpp += relative-gen relativecm-gen
module_programs_cpp_test := relative_rme_test relativecm_rme_test
module_programs_cpp_test += wigner_test
# module_programs_f :=
# module_generated :=

################################################################
# library creation flag
################################################################

$(eval $(library))

$(eval $(end-module))
