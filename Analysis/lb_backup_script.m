%lb_backup_script

% Clark_P30 (LB and ATN projects)
backup({'F:\ClarkP30_Recordings'},'Z:\Laura_Berkowitz\Projects')

%% ATNppc
% push to server
backup({'F:\Projects\ATNppc'},'Z:\ryan_harvey')
% pull from server
backup({'Z:\ryan_harvey\Projects\ATNppc'},'F:','ignore_path','ryan_harvey')