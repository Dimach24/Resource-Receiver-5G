classdef ResourceReceiver
    methods(Static)
        function [pbch,dmrs]=PbchExtraction(Rgrid,toffset,foffset,NCellId)
            % extracts pbch bitstream and dmrs complex 
            % amplitudes from the resource grid
            
            nu=mod(NCellId,4);
            
            % subcarrier index initialization 
            d_solid_i=(0:4:236)+nu;
            p_solid_i=0:1:239;
            p_solid_i(d_solid_i+1)=[];

            d_splitted_i=[(0:4:44)+nu, (192:4:236)+nu];
            p_splitted_i=[0:47,192:239];
            p_splitted_i((0:4:92)+nu+1)=[];

            % extraction from the resource grid 
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
            % converting into bitstream
            pbch=QpskDemodulation(pbch);
        end
        function blockIndexLsb=PbchDmRsProcessing(dmrs_linearized,NCellId)
            % finds the block index that corresponds to the reference signal

            dmrs_bank=zeros(8,144);
            % generating signals
            for i=0:7
                dmrs_bank(i+1,:)=generatePbchDmRs(i,NCellId);
            end
            % correlating signals
            for i=1:8
                corr_data(i,:)=abs(xcorr(dmrs_bank(i,:),dmrs_linearized));
            end
            % maximum likehood search
            corr_max=max(corr_data,[],2);
            [~,blockIndexLsb]=max(corr_max);
            blockIndexLsb=blockIndexLsb-1;
        end
        function out_seq = Descramble(in_seq, N_Cell_ID, L_max, block_index)
            %ScrambleProcedure of revererse scrambling
            % after demodulation [7.3.3.1, TS 38.211]
            arguments
                in_seq (1,:) % input sequence (boolean matrix)
                N_Cell_ID (1,1)
                L_max (1,1) % maximum number of candidate SS/PBCH blocks in half frame [5, TS 38.213]
                block_index (1,1) % candidate SS/PBCH block index
            end
            
            %init
            A = length(in_seq);
            s = zeros(1,A);
            M = A;
            
            %determinaton of nu
            block_index = fliplr(dec2bin(block_index,3));
            if L_max == 4
                nu = [block_index(2) block_index(1)];
            else
                nu = [block_index(3) block_index(2) block_index(1)];
            end
            nu = bin2dec(num2str(nu));
            
            %determination of c
            x1 = zeros(1,2000);
            x2 = zeros(1,2000);
            x1(1) = 1;
            x1(2:31) = 0;
            x2(1:31) = fliplr(int2bit(N_Cell_ID,31)); %c_init = N_Cell_ID
            for n = 1:2000
                x1(n+31) = mod(x1(n+3)+x1(n),2);
                x2(n+31) = mod(x2(n+3)+x2(n+2)+x2(n+1)+x2(n),2);
                n1 = 1:160;
                c(n1) = mod(x1(n1+1600)+x2(n1+1600),2);
            end
            
            %determination of s
            i = 0;
            j = 0;
            while i < A
                s(1+i) = c(1+mod(j+nu*M,160));
                j = j+1;
                i = i+1;
            end
            
            %descrambling procedure
            out_seq = zeros (1,A);
            for i = 1:A
                out_seq(i) = mod(in_seq(i)+ s(i),2);
            end
        end
        function [bitstream,i_ssb_lsb]=getBitstream(Rgrid,toffset,foffset,NCellId,L_max)
            [pbch,dmrs]=ResourceReceiver.PbchExtraction(Rgrid,toffset,foffset,NCellId);
            i_ssb_lsb=ResourceReceiver.PbchDmRsProcessing(dmrs,NCellId);
            bitstream=ResourceReceiver.Descramble(pbch,NCellId,L_max,i_ssb_lsb);
        end
    end
end