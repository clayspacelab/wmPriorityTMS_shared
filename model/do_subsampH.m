function [sub_data] = do_subsampH(data,subjid,permidx,condVec,maxperms)
%helper function for subsampling high data from a fixed (per subject)
%permutation list to match amount of data in low priority cond.

if ~exist('maxperms','var')
    maxperms = 20;
end
if ~exist('condVec','var')
    condVec = fieldnames(data);
end

%set random seed
s = RandStream('mt19937ar','Seed',subjid);


%get data lengths
%high_len = structfun(@(x) length(x{1}),data);
%low_len = structfun(@(x) length(x{2}),data);

perms = struct();
for i=1:length(condVec)
    %don't count NaNs toward totals
    high_len = sum(~isnan(data.(condVec{i}){1}));
    low_len = sum(~isnan(data.(condVec{i}){2}));

    perms.(condVec{i}) = zeros(low_len,maxperms);
    for j=1:maxperms
        perms.(condVec{i})(:,j) = randperm(s,high_len,low_len);
    end
end

%note that we are not going to touch data not in condVec, but we're not going to remove it either 
sub_data = data;

for f=1:length(condVec)
    %remove NaNs
    this_data = sub_data.(condVec{f}){1};
    this_data(isnan(this_data)) = [];
    %save perm subset
    sub_data.(condVec{f}){1} = this_data(perms.(condVec{f})(:,permidx));
end


end