function bitstream=QpskDemodulation(sequence)
    bitstream=zeros(1,length(sequence)*2,"uint8");
    for i=1:length(sequence)
        bitstream(2*i-1)=(1-sign(real(sequence(i))))/2;
        bitstream(2*i)=(1-sign(imag(sequence(i))))/2;
    end
end