function [slt, W] = sero_seq_seq2sliceThickness(seq, dz)
% function [slt, W] = sero_seq_seq2sliceThickness(seq, dz)

nEvents = numel(seq.blockEvents);

Lib = {'initialized'};

c = 1;
for i = 1:nEvents

    cEvent = seq.blockEvents{i};

    if cEvent(2) && cEvent(5) && cEvent(7)

        block = seq.getBlock(i);

        rfWF = block.rf.signal';
        rfT  = block.rf.t'*1000;
        gzA  = mr.convert(block.gz.waveform(2), 'Hz/m', 'mT/m');

        sha = sero_seq_array2sha([abs(rfWF(:)); angle(rfWF(:)); rfT(:); gzA]);

        indLib = find(ismember(Lib(:,1),{sha}));

        if isempty(indLib)
            [M, z, sliceThick] = sero_seq_sliceProfile(rfWF, rfT, 90, gzA, [], 1);

            Mc  = M(1,:) + M(2,:)*1i;

            weight = sero_recon_sliceProf2weight(Mc, z, dz, sliceThick, 0);

            indEnd = size(Lib,1)+1;

            Lib{indEnd, 1} = sha;
            Lib{indEnd, 2} = sliceThick;
            Lib{indEnd, 3} = weight;

        else
            sliceThick = Lib{indLib, 2};
            weight = Lib{indLib, 3};

        end

        slt(c) = sliceThick;
        W{c}   = weight;

        disp(c);
        c = c+1;        
    end

end
