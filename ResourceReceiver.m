classdef ResourceReceiver
    methods(Static)
        function [pbch,dmrs]=PbchExtraction(Rgrid,toffset,foffset,NCellId)
            nu=mod(NCellId,4);
            indexes_fsolid=find(mod(1:1:240,4)==1);
            indexes_fsplitted=indexes_fsolid(indexes_fsolid<46 | indexes_fsolid>192);
            indexes_fsolid=indexes_fsoled + nu +foffset;
            indexes_fsplitted=indexes_fsplitted + nu +foffset;
            dmrs=[...
                reshape(Rgrid(indexes_fsolid,   toffset+1),[],1),...
                reshape(Rgrid(indexes_fsplitted,toffset+2),[],1),...
                reshape(Rgrid(indexes_fsolid,   toffset+3),[],1)...
                ];
            pbch=[...
                reshape(Rgrid(setdiff(1:240+foffset,indexes_fsolid)     ,toffset+1),[],1),...
                reshape(Rgrid(setdiff(1:48+foffset,indexes_fsplitted)   ,toffset+2),[],1),...
                reshape(Rgrid(setdiff(193:240+foffset,indexes_fsplitted),toffset+2),[],1),...
                reshape(Rgrid(setdiff(1:240+foffset,indexes_fsolid)     ,toffset+3),[],1),...
                ];
            pbch=QpskDemodulation(pbch);
        end
        function blockIndexLsb=PbchDmRsProcessing(dmrs_linearized,NCellId)
            for i=0:7
                dmrs_bank(i+1,:)=generatePbchDmRs(i,NCellId);
            end
            corr_data=xcorr(dmrs_bank,dmrs_linearized);
            corr_max=max(corr_data);
            [~,blockIndexLsb]=max(corr_max,[]);
        end
    end
end