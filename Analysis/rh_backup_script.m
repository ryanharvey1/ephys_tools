%rh_backup_script



%% HPCatn
backup({'F:\Projects\HPCatn'},'Z:\ryan_harvey')



%% PAE
backup({'D:\Projects\PAE_PlaceCell'},'Z:\ryan_harvey')



%% ATNppc
% push to server
backup({'F:\Projects\ATNppc'},'Z:\ryan_harvey')
% pull from server
backup({'Z:\ryan_harvey\Projects\ATNppc'},'F:','ignore_path','ryan_harvey')


