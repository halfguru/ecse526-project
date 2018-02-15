function create_parpool(clusterName, nWorkers, forceRestart)

% clusterName = 'local';
% nWorkers = 7;
poolobj = gcp('nocreate');
if isempty(poolobj)
    emptyPool = true;
    numWorkers = 0;
else
    emptyPool = false;
    numWorkers = poolobj.NumWorkers;
end
if (emptyPool) || (numWorkers ~= nWorkers) || (forceRestart)
    delete(poolobj);
    cc = parcluster(clusterName);
    if strcmpi(clusterName,'local')
        cc.NumWorkers = nWorkers;
        parpool(cc,cc.NumWorkers);
    else
        parpool(cc);
    end
end

end