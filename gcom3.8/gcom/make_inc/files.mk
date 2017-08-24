#                        GCOM Makefile Include File

#                   This file should not need to be edited

# *****************************COPYRIGHT*******************************
# (c) CROWN COPYRIGHT 2001, Met Office, All Rights Reserved. 
# Please refer to Copyright file in top level GCOM directory
#                 for further details
# *****************************COPYRIGHT*******************************

GC_SRC= gcom_mod.F90 \
gc_globals_mod.F90 \
gc__errlim.F90 \
gc__get_mpi_type.F90 \
gc__init_internals.F90 \
gc__init_timer.F90 \
gc__stamp.F90 \
gc__timer_funcs.F90 \
gc_abort.F90 \
gc_bbcast.F90 \
gc_brecv.F90 \
gc_bsend.F90 \
gc_cbcast.F90 \
gc_config.F90 \
gc_crecv.F90 \
gc_csend.F90 \
gc_exit.F90 \
gc_get_communicator.F90 \
gc_getopt.F90 \
gc_gsync.F90 \
gc_ibcast.F90 \
gc_imax.F90 \
gc_imin.F90 \
gc_init.F90 \
gc_irecv.F90 \
gc_isend.F90 \
gc_isum.F90 \
gc_me.F90 \
gc_nproc.F90 \
gc_rbcast.F90 \
gc_rmax.F90 \
gc_rmin.F90 \
gc_rrecv.F90 \
gc_rsend.F90 \
gc_rsum.F90 \
gc_rsumr.F90 \
gc_set_communicator.F90 \
gc_setopt.F90 \
gc_ssync.F90 \
gc__flush.F90


GC_COMMON_INC= \
        ${INCLUDE_DIR}/gc__init_common.h \
        ${INCLUDE_DIR}/gc__mpi_common.h \
        ${INCLUDE_DIR}/gc__mpi_types_common.h

GC_PROLOG_INC= \
${INCLUDE_DIR}/gc_aix.h \
${INCLUDE_DIR}/gc_assorted.h \
${INCLUDE_DIR}/gc_com.h \
${INCLUDE_DIR}/gc_constants.h \
${INCLUDE_DIR}/gc_end_timer.h \
${INCLUDE_DIR}/gc_functions.h \
${INCLUDE_DIR}/gc_interface.h \
${INCLUDE_DIR}/gc_kinds.h \
${INCLUDE_DIR}/gc_limits.h \
${INCLUDE_DIR}/gc_mtags.h \
${INCLUDE_DIR}/gc_precision.h \
${INCLUDE_DIR}/gc_prolog.h \
${INCLUDE_DIR}/gc_set_precision.h \
${INCLUDE_DIR}/gc_start_timer.h \
${INCLUDE_DIR}/gc_timer.h

GCG_SRC= gcg__errlim.F90 \
         gcg__init_internals.F90 \
         gcg__mpi_rank.F90 \
         gcg_config.F90 \
         gcg_ibcast.F90 \
         gcg_imax.F90 \
         gcg_imin.F90 \
         gcg_isum.F90 \
         gcg_me.F90 \
         gcg_ralltoalle.F90 \
         gcg_rbcast.F90 \
         gcg_rmax.F90 \
         gcg_rmin.F90 \
         gcg_rsum.F90 \
         gcg_rsumr.F90 \
         gcg_rvecshift.F90 \
         gcg_rvecsumf.F90 \
         gcg_rvecsumr.F90 \
         gcg_split.F90 \
         gcg_ssync.F90 \
         gcg_ralltoalle_multi.F90 \
         gcg_rvecsumrf.F90

GCG_COMMON_INC= \
         ${INCLUDE_DIR}/gcg__barrier_common.h \
         ${INCLUDE_DIR}/gcg__grmpi_common.h \
         ${INCLUDE_DIR}/gcg__grstore_common.h \
         ${INCLUDE_DIR}/gcg__rotate_common.h \
         ${INCLUDE_DIR}/gcg__split_common.h 

GCG_PROLOG_INC= \
          ${GC_PROLOG_INC} \
          ${INCLUDE_DIR}/gcg_constants.h \
          ${INCLUDE_DIR}/gcg_limits.h \
          ${INCLUDE_DIR}/gcg_mtags.h \
          ${INCLUDE_DIR}/gcg_prolog.h
          

MPL_SRC= mpl.F90 \
mpl_abort.F90 \
mpl_allgather.F90 \
mpl_allgatherv.F90 \
mpl_allreduce.F90 \
mpl_alltoall.F90 \
mpl_alltoallv.F90 \
mpl_barrier.F90 \
mpl_bcast.F90 \
mpl_bsend.F90 \
mpl_buffer_attach.F90 \
mpl_comm_call_errhandler.F90 \
mpl_comm_create.F90 \
mpl_comm_create_errhandler.F90 \
mpl_comm_dup.F90 \
mpl_comm_get_attr.F90 \
mpl_comm_get_errhandler.F90 \
mpl_comm_group.F90 \
mpl_comm_rank.F90 \
mpl_comm_set_errhandler.F90 \
mpl_comm_size.F90 \
mpl_comm_split.F90 \
mpl_error_class.F90 \
mpl_error_string.F90 \
mpl_finalize.F90 \
mpl_finalized.F90 \
mpl_gather.F90 \
mpl_gatherv.F90 \
mpl_get.F90 \
mpl_group_free.F90 \
mpl_group_incl.F90 \
mpl_group_translate_ranks.F90 \
mpl_init.F90 \
mpl_init_thread.F90 \
mpl_initialized.F90 \
mpl_iprobe.F90 \
mpl_irecv.F90 \
mpl_isend.F90 \
mpl_probe.F90 \
mpl_query_thread.F90 \
mpl_recv.F90 \
mpl_reduce.F90 \
mpl_scatterv.F90 \
mpl_send.F90 \
mpl_sendrecv.F90 \
mpl_ssend.F90 \
mpl_test.F90 \
mpl_testall.F90 \
mpl_type_commit.F90 \
mpl_type_create_resized.F90 \
mpl_type_extent.F90 \
mpl_type_free.F90 \
mpl_type_struct.F90 \
mpl_type_vector.F90 \
mpl_wait.F90 \
mpl_waitall.F90 \
mpl_win_create.F90 \
mpl_win_fence.F90 \
mpl_win_free.F90 
