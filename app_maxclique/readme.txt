1. mc_run.cpp: original mcf program, with sequential 1st round split. no time delay divide conquer and no expand split. VNUM_ALLOWED_BEFORE_SPLIT is input argument
2. mc_tddq_run.cpp: mcf program with sequential 1st round split and time delay divide conquer. no expand split
3. tddq_split_run.cpp: mcf program with sequential 1st round split and time delay divide conquer and expand split. VNUM_ALLOWED_BEFORE_SPLIT 400000
4. tddq_seq_split_run.cpp: mcf program with sequential 1st round split and time delay divide conquer and expand split. VNUM_ALLOWED_BEFORE_SPLIT is input argument
5. tddq_para_split_run.cpp: mcf program with parallel 1st round split and time delay divide conquer and expand split. VNUM_ALLOWED_BEFORE_SPLIT is input argument
6. tddq_non_expand_run.cpp: mcf program with sequential 1st round split and time delay divide conquer. No expand split. VNUM_ALLOWED_BEFORE_SPLIT is input argument
7. tddq_split_run_only.cpp: mcf program with NO 1st round split. There is time delay divide conquer and expand split. 