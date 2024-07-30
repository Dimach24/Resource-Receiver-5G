classdef ResourceReceiver
    methods(Static)
        function [pbch,dmrs]=PbchExtraction(Rgrid,toffset,foffset,NCellId)
            nu=mod(NCellId,4);
            d_solid_i=(0:4:236)+nu;
            p_solid_i=0:1:239;
            p_solid_i(d_solid_i+1)=[];

            d_splitted_i=[(0:4:44)+nu, (192:4:236)+nu];
            p_splitted_i=[0:47,192:239];
            p_splitted_i((0:4:92)+nu+1)=[];

            
            dmrs=[...
                Rgrid(foffset+d_solid_i+1,toffset+1+1).',...
                Rgrid(foffset+d_splitted_i+1,toffset+2+1).',...
                Rgrid(foffset+d_solid_i+1,toffset+3+1).'...
                ];
            pbch=[...
                Rgrid(foffset+p_solid_i+1,toffset+1+1).',...
                Rgrid(foffset+p_splitted_i+1,toffset+2+1).',...
                Rgrid(foffset+p_solid_i+1,toffset+3+1).'...
                ];
            % Rgrid(foffset+p_solid_i+1,toffset+1+1)=ones(1,180)*10;
            % Rgrid(foffset+p_splitted_i+1,toffset+2+1)=ones(1,72)*10;
            % Rgrid(foffset+p_solid_i+1,toffset+3+1)=ones(1,180)*10;
            % Rgrid(foffset+d_solid_i+1,toffset+1+1)=ones(1,60)*5;
            % Rgrid(foffset+d_splitted_i+1,toffset+2+1)=ones(1,24)*5;
            % Rgrid(foffset+d_solid_i+1,toffset+3+1)=ones(1,60)*5;
            % plt=pcolor(abs(Rgrid(1:301,1:end)));
            % plt.EdgeColor='none';
            % ca=gca();
            % ca.YDir='normal';
            % xlim([1,50]);
            % xlabel('l+1 (номер OFDM символа +1)')
            % ylabel('k (номер поднесущей)')

            pbch=QpskDemodulation(pbch);
        end
        function blockIndexLsb=PbchDmRsProcessing(dmrs_linearized,NCellId)
            for i=0:7
                dmrs_bank(i+1,:)=generatePbchDmRs(i,NCellId);
            end
            for i=1:8
                corr_data(i,:)=abs(xcorr(dmrs_bank(i,:),dmrs_linearized));
            end
            corr_max=max(corr_data,[],2);
            [~,blockIndexLsb]=max(corr_max);
            blockIndexLsb=blockIndexLsb-1;
        end
    end
end