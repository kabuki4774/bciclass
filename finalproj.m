clear all

%load all data collected
proj = load('bciprojectdata.mat');

trials = 32;
channels = 16;
fs = 1200;

%simple filter for later
fc = 2;
[b,a] = butter(4,fc/(fs/2),'high');
freqz(b,a);
figure;

%average observation for noise reduction

%bahram - this is unused as the protocol for this data has been unclear.

for ch = 1:channels
    %collect data from all trials
    for tr = 1:trials
        proj.bahram.consec.cumm.temp(:,tr) = proj.bahram.consec.epoch_data_new_zero1{1,tr}(:,ch);
    end
    %average trials from all data in one channel
    proj.bahram.consec.avg1(:,ch) = mean(proj.bahram.consec.cumm.temp,2);
    clear proj.bahram.consec.cumm.temp;
end

for ch = 1:channels
    %collect data from all trials
    for tr = 1:trials
        proj.bahram.consec.cumm.temp(:,tr) = proj.bahram.consec.epoch_data_new_zero2{1,tr}(:,ch);
    end
    %average trials from all data in one channel
    proj.bahram.consec.avg2(:,ch) = mean(proj.bahram.consec.cumm.temp,2);
    clear proj.bahram.consec.cumm.temp;
end

for ch = 1:channels
    %collect data from all trials
    for tr = 1:trials
        proj.bahram.consec.cumm.temp(:,tr) = proj.bahram.consec.epoch_data_new_zero4{1,tr}(:,ch);
    end
    %average trials from all data in one channel
    proj.bahram.consec.avg4(:,ch) = mean(proj.bahram.consec.cumm.temp,2);
    clear proj.bahram.consec.cumm.temp;
end

%sara

for ch = 1:channels
    %collect data from all trials
    for tr = 1:trials
        proj.sara.consec.cumm.temp(:,tr) = proj.sara.consec.epoch_data_new_zero1{1,tr}(:,ch);
    end
    %average trials from all data in one channel
    proj.sara.consec.avg1(:,ch) = mean(proj.sara.consec.cumm.temp,2);
    clear proj.sara.consec.cumm.temp;
end

for ch = 1:channels
    %collect data from all trials
    for tr = 1:trials
        proj.sara.consec.cumm.temp(:,tr) = proj.sara.consec.epoch_data_new_zero2{1,tr}(:,ch);
    end
    %average trials from all data in one channel
    proj.sara.consec.avg2(:,ch) = mean(proj.sara.consec.cumm.temp,2);
    clear proj.sara.consec.cumm.temp;
end

for ch = 1:channels
    %collect data from all trials
    for tr = 1:trials
        proj.sara.consec.cumm.temp(:,tr) = proj.sara.consec.epoch_data_new_zero4{1,tr}(:,ch);
    end
    %average trials from all data in one channel
    proj.sara.consec.avg4(:,ch) = mean(proj.sara.consec.cumm.temp,2);
    clear proj.sara.consec.cumm.temp;
end

for ch = 1:channels
    %collect data from all trials
    for tr = 1:trials
        proj.sara.sim.cumm.temp(:,tr) = proj.sara.sim.epoch_data_new_zero1{1,tr}(:,ch);
    end
    %average trials from all data in one channel
    proj.sara.sim.avg1(:,ch) = mean(proj.sara.sim.cumm.temp,2);
    clear proj.sara.sim.cumm.temp;
end

for ch = 1:channels
    %collect data from all trials
    for tr = 1:trials
        proj.sara.sim.cumm.temp(:,tr) = proj.sara.sim.epoch_data_new_zero2{1,tr}(:,ch);
    end
    %average trials from all data in one channel
    proj.sara.sim.avg2(:,ch) = mean(proj.sara.sim.cumm.temp,2);
    clear proj.sara.sim.cumm.temp;
end

for ch = 1:channels
    %collect data from all trials
    for tr = 1:trials
        proj.sara.sim.cumm.temp(:,tr) = proj.sara.sim.epoch_data_new_zero4{1,tr}(:,ch);
    end
    %average trials from all data in one channel
    proj.sara.sim.avg4(:,ch) = mean(proj.sara.sim.cumm.temp,2);
    clear proj.sara.sim.cumm.temp;
end

%roohi

for ch = 1:channels
    %collect data from all trials
    for tr = 1:trials
        proj.roohi.consec.cumm.temp(:,tr) = proj.roohi.consec.epoch_data_new_zero1{1,tr}(:,ch);
    end
    %average trials from all data in one channel
    proj.roohi.consec.avg1(:,ch) = mean(proj.roohi.consec.cumm.temp,2);
    clear proj.roohi.consec.cumm.temp;
end

for ch = 1:channels
    %collect data from all trials
    for tr = 1:trials
        proj.roohi.consec.cumm.temp(:,tr) = proj.roohi.consec.epoch_data_new_zero2{1,tr}(:,ch);
    end
    %average trials from all data in one channel
    proj.roohi.consec.avg2(:,ch) = mean(proj.roohi.consec.cumm.temp,2);
    clear proj.roohi.consec.cumm.temp;
end

for ch = 1:channels
    %collect data from all trials
    for tr = 1:trials
        proj.roohi.consec.cumm.temp(:,tr) = proj.roohi.consec.epoch_data_new_zero4{1,tr}(:,ch);
    end
    %average trials from all data in one channel
    proj.roohi.consec.avg4(:,ch) = mean(proj.roohi.consec.cumm.temp,2);
    clear proj.roohi.consec.cumm.temp;
end

for ch = 1:channels
    %collect data from all trials
    for tr = 1:trials
        proj.roohi.sim.cumm.temp(:,tr) = proj.roohi.sim.epoch_data_new_zero1{1,tr}(:,ch);
    end
    %average trials from all data in one channel
    proj.roohi.sim.avg1(:,ch) = mean(proj.roohi.sim.cumm.temp,2);
    clear proj.roohi.sim.cumm.temp;
end

for ch = 1:channels
    %collect data from all trials
    for tr = 1:trials
        proj.roohi.sim.cumm.temp(:,tr) = proj.roohi.sim.epoch_data_new_zero2{1,tr}(:,ch);
    end
    %average trials from all data in one channel
    proj.roohi.sim.avg2(:,ch) = mean(proj.roohi.sim.cumm.temp,2);
    clear proj.roohi.sim.cumm.temp;
end

for ch = 1:channels
    %collect data from all trials
    for tr = 1:trials
        proj.roohi.sim.cumm.temp(:,tr) = proj.roohi.sim.epoch_data_new_zero4{1,tr}(:,ch);
    end
    %average trials from all data in one channel
    proj.roohi.sim.avg4(:,ch) = mean(proj.roohi.sim.cumm.temp,2);
    clear proj.roohi.sim.cumm.temp;
end

%tanya

for ch = 1:channels
    %collect data from all trials
    for tr = 1:trials
        proj.tanya.consec.cumm.temp(:,tr) = proj.tanya.consec.epoch_data_new_zero1{1,tr}(:,ch);
    end
    %average trials from all data in one channel
    proj.tanya.consec.avg1(:,ch) = mean(proj.tanya.consec.cumm.temp,2);
    clear proj.tanya.consec.cumm.temp;
end

for ch = 1:channels
    %collect data from all trials
    for tr = 1:trials
        proj.tanya.consec.cumm.temp(:,tr) = proj.tanya.consec.epoch_data_new_zero2{1,tr}(:,ch);
    end
    %average trials from all data in one channel
    proj.tanya.consec.avg2(:,ch) = mean(proj.tanya.consec.cumm.temp,2);
    clear proj.tanya.consec.cumm.temp;
end

for ch = 1:channels
    %collect data from all trials
    for tr = 1:trials
        proj.tanya.consec.cumm.temp(:,tr) = proj.tanya.consec.epoch_data_new_zero4{1,tr}(:,ch);
    end
    %average trials from all data in one channel
    proj.tanya.consec.avg4(:,ch) = mean(proj.tanya.consec.cumm.temp,2);
    clear proj.tanya.consec.cumm.temp;
end

for ch = 1:channels
    %collect data from all trials
    for tr = 1:trials
        proj.tanya.sim.cumm.temp(:,tr) = proj.tanya.sim.epoch_data_new_zero1{1,tr}(:,ch);
    end
    %average trials from all data in one channel
    proj.tanya.sim.avg1(:,ch) = mean(proj.tanya.sim.cumm.temp,2);
    clear proj.tanya.sim.cumm.temp;
end

for ch = 1:channels
    %collect data from all trials
    for tr = 1:trials
        proj.tanya.sim.cumm.temp(:,tr) = proj.tanya.sim.epoch_data_new_zero2{1,tr}(:,ch);
    end
    %average trials from all data in one channel
    proj.tanya.sim.avg2(:,ch) = mean(proj.tanya.sim.cumm.temp,2);
    clear proj.tanya.sim.cumm.temp;
end

for ch = 1:channels
    %collect data from all trials
    for tr = 1:trials
        proj.tanya.sim.cumm.temp(:,tr) = proj.tanya.sim.epoch_data_new_zero4{1,tr}(:,ch);
    end
    %average trials from all data in one channel
    proj.tanya.sim.avg4(:,ch) = mean(proj.tanya.sim.cumm.temp,2);
    clear proj.tanya.sim.cumm.temp;
end

%take averages and perform cleaning
%   - filter
%   - whiten
%   - normalize
%take psd for visualizing vibration of motors

for ch = 1:channels
    proj.tanya.consec.avg1(:,ch)=filter(b,a,proj.tanya.consec.avg1(:,ch));
    %whiten data
    fudgefactor = 0.0001;
    proj.tanya.consec.avg1(:,ch) = bsxfun(@minus, proj.tanya.consec.avg1(:,ch), mean(proj.tanya.consec.avg1(:,ch)));
    A = proj.tanya.consec.avg1(:,ch)'*proj.tanya.consec.avg1(:,ch);
    [V,D] = eig(A);
    proj.tanya.consec.avg1(:,ch) = proj.tanya.consec.avg1(:,ch)*V*diag(1./(diag(D)+fudgefactor).^(1/2))*V';
    %normalize
    proj.tanya.consec.avg1(:,ch)=proj.tanya.consec.avg1(:,ch)-mean(proj.tanya.consec.avg1(:,ch));
    proj.tanya.consec.avg1(:,ch)=proj.tanya.consec.avg1(:,ch)/std(proj.tanya.consec.avg1(:,ch));
    %find PSD
    [ob2.tanya.ps2(:,ch),ob2.tanya.freq2(:,ch)] = pwelch(proj.tanya.consec.avg1(:,ch),fs,fs/2);
    subplot(8,2,ch)
    title(num2str(ch))
    plot(ob2.tanya.freq2(:,ch)/(2*pi)*fs,ob2.tanya.ps2(:,ch))
    xlim([0 40])
    % computing of the peak (crest) factor
    ob2.tanya.rms2(ch) = sqrt(mean(ob2.tanya.ps2(:,ch).^2));
    ob2.tanya.peak2(ch) = max(abs(ob2.tanya.ps2(:,ch)));
    ob2.tanya.Q2(ch) = 20*log10(ob2.tanya.peak2(ch)/ob2.tanya.rms2(ch));
    disp(['Peak (crest) factor Q = ' num2str(ob2.tanya.Q2(ch)) ' dB'])
    %spectrogram(ps1(:,ch), [], [], [], fs, 'yaxis')
end

for ch = 1:channels
    proj.roohi.consec.avg1(:,ch)=filter(b,a,proj.roohi.consec.avg1(:,ch));
    %whiten data
    fudgefactor = 0.0001;
    proj.roohi.consec.avg1(:,ch) = bsxfun(@minus, proj.roohi.consec.avg1(:,ch), mean(proj.roohi.consec.avg1(:,ch)));
    A = proj.roohi.consec.avg1(:,ch)'*proj.roohi.consec.avg1(:,ch);
    [V,D] = eig(A);
    proj.roohi.consec.avg1(:,ch) = proj.roohi.consec.avg1(:,ch)*V*diag(1./(diag(D)+fudgefactor).^(1/2))*V';
    %normalize
    proj.roohi.consec.avg1(:,ch)=proj.roohi.consec.avg1(:,ch)-mean(proj.roohi.consec.avg1(:,ch));
    proj.roohi.consec.avg1(:,ch)=proj.roohi.consec.avg1(:,ch)/std(proj.roohi.consec.avg1(:,ch));
    %find PSD
    [ob2.roohi.ps2(:,ch),ob2.roohi.freq2(:,ch)] = pwelch(proj.roohi.consec.avg1(:,ch),fs,fs/2);
    subplot(8,2,ch)
    title(num2str(ch))
    plot(ob2.roohi.freq2(:,ch)/(2*pi)*fs,ob2.roohi.ps2(:,ch))
    xlim([0 40])
    % computing of the peak (crest) factor
    ob2.roohi.rms2(ch) = sqrt(mean(ob2.roohi.ps2(:,ch).^2));
    ob2.roohi.peak2(ch) = max(abs(ob2.roohi.ps2(:,ch)));
    ob2.roohi.Q2(ch) = 20*log10(ob2.roohi.peak2(ch)/ob2.roohi.rms2(ch));
    disp(['Peak (crest) factor Q = ' num2str(ob2.roohi.Q2(ch)) ' dB'])
    %spectrogram(ps1(:,ch), [], [], [], fs, 'yaxis')
end

for ch = 1:channels
    proj.sara.consec.avg1(:,ch)=filter(b,a,proj.sara.consec.avg1(:,ch));
    %whiten data
    fudgefactor = 0.0001;
    proj.sara.consec.avg1(:,ch) = bsxfun(@minus, proj.sara.consec.avg1(:,ch), mean(proj.sara.consec.avg1(:,ch)));
    A = proj.sara.consec.avg1(:,ch)'*proj.sara.consec.avg1(:,ch);
    [V,D] = eig(A);
    proj.sara.consec.avg1(:,ch) = proj.sara.consec.avg1(:,ch)*V*diag(1./(diag(D)+fudgefactor).^(1/2))*V';
    %normalize
    proj.sara.consec.avg1(:,ch)=proj.sara.consec.avg1(:,ch)-mean(proj.sara.consec.avg1(:,ch));
    proj.sara.consec.avg1(:,ch)=proj.sara.consec.avg1(:,ch)/std(proj.sara.consec.avg1(:,ch));
    %find PSD
    [ob2.sara.ps2(:,ch),ob2.sara.freq2(:,ch)] = pwelch(proj.sara.consec.avg1(:,ch),fs,fs/2);
    subplot(8,2,ch)
    title(num2str(ch))
    plot(ob2.sara.freq2(:,ch)/(2*pi)*fs,ob2.sara.ps2(:,ch))
    xlim([0 40])
    % computing of the peak (crest) factor
    ob2.sara.rms2(ch) = sqrt(mean(ob2.sara.ps2(:,ch).^2));
    ob2.sara.peak2(ch) = max(abs(ob2.sara.ps2(:,ch)));
    ob2.sara.Q2(ch) = 20*log10(ob2.sara.peak2(ch)/ob2.sara.rms2(ch));
    disp(['Peak (crest) factor Q = ' num2str(ob2.sara.Q2(ch)) ' dB'])
    %spectrogram(ps1(:,ch), [], [], [], fs, 'yaxis')
end

for ch = 1:channels
    proj.tanya.sim.avg1(:,ch)=filter(b,a,proj.tanya.sim.avg1(:,ch));
    %whiten data
    fudgefactor = 0.0001;
    proj.tanya.sim.avg1(:,ch) = bsxfun(@minus, proj.tanya.sim.avg1(:,ch), mean(proj.tanya.sim.avg1(:,ch)));
    A = proj.tanya.sim.avg1(:,ch)'*proj.tanya.sim.avg1(:,ch);
    [V,D] = eig(A);
    proj.tanya.sim.avg1(:,ch) = proj.tanya.sim.avg1(:,ch)*V*diag(1./(diag(D)+fudgefactor).^(1/2))*V';
    %normalize
    proj.tanya.sim.avg1(:,ch)=proj.tanya.sim.avg1(:,ch)-mean(proj.tanya.sim.avg1(:,ch));
    proj.tanya.sim.avg1(:,ch)=proj.tanya.sim.avg1(:,ch)/std(proj.tanya.sim.avg1(:,ch));
    %find PSD
    [ob2.tanya.ps2(:,ch),ob2.tanya.freq2(:,ch)] = pwelch(proj.tanya.sim.avg1(:,ch),fs,fs/2);
    subplot(8,2,ch)
    title(num2str(ch))
    plot(ob2.tanya.freq2(:,ch)/(2*pi)*fs,ob2.tanya.ps2(:,ch))
    xlim([0 40])
    % computing of the peak (crest) factor
    ob2.tanya.rms2(ch) = sqrt(mean(ob2.tanya.ps2(:,ch).^2));
    ob2.tanya.peak2(ch) = max(abs(ob2.tanya.ps2(:,ch)));
    ob2.tanya.Q2(ch) = 20*log10(ob2.tanya.peak2(ch)/ob2.tanya.rms2(ch));
    disp(['Peak (crest) factor Q = ' num2str(ob2.tanya.Q2(ch)) ' dB'])
    %spectrogram(ps1(:,ch), [], [], [], fs, 'yaxis')
end

for ch = 1:channels
    proj.roohi.sim.avg1(:,ch)=filter(b,a,proj.roohi.sim.avg1(:,ch));
    %whiten data
    fudgefactor = 0.0001;
    proj.roohi.sim.avg1(:,ch) = bsxfun(@minus, proj.roohi.sim.avg1(:,ch), mean(proj.roohi.sim.avg1(:,ch)));
    A = proj.roohi.sim.avg1(:,ch)'*proj.roohi.sim.avg1(:,ch);
    [V,D] = eig(A);
    proj.roohi.sim.avg1(:,ch) = proj.roohi.sim.avg1(:,ch)*V*diag(1./(diag(D)+fudgefactor).^(1/2))*V';
    %normalize
    proj.roohi.sim.avg1(:,ch)=proj.roohi.sim.avg1(:,ch)-mean(proj.roohi.sim.avg1(:,ch));
    proj.roohi.sim.avg1(:,ch)=proj.roohi.sim.avg1(:,ch)/std(proj.roohi.sim.avg1(:,ch));
    %find PSD
    [ob2.roohi.ps2(:,ch),ob2.roohi.freq2(:,ch)] = pwelch(proj.roohi.sim.avg1(:,ch),fs,fs/2);
    subplot(8,2,ch)
    title(num2str(ch))
    plot(ob2.roohi.freq2(:,ch)/(2*pi)*fs,ob2.roohi.ps2(:,ch))
    xlim([0 40])
    % computing of the peak (crest) factor
    ob2.roohi.rms2(ch) = sqrt(mean(ob2.roohi.ps2(:,ch).^2));
    ob2.roohi.peak2(ch) = max(abs(ob2.roohi.ps2(:,ch)));
    ob2.roohi.Q2(ch) = 20*log10(ob2.roohi.peak2(ch)/ob2.roohi.rms2(ch));
    disp(['Peak (crest) factor Q = ' num2str(ob2.roohi.Q2(ch)) ' dB'])
    %spectrogram(ps1(:,ch), [], [], [], fs, 'yaxis')
end

for ch = 1:channels
    proj.sara.sim.avg1(:,ch)=filter(b,a,proj.sara.sim.avg1(:,ch));
    %whiten data
    fudgefactor = 0.0001;
    proj.sara.sim.avg1(:,ch) = bsxfun(@minus, proj.sara.sim.avg1(:,ch), mean(proj.sara.sim.avg1(:,ch)));
    A = proj.sara.sim.avg1(:,ch)'*proj.sara.sim.avg1(:,ch);
    [V,D] = eig(A);
    proj.sara.sim.avg1(:,ch) = proj.sara.sim.avg1(:,ch)*V*diag(1./(diag(D)+fudgefactor).^(1/2))*V';
    %normalize
    proj.sara.sim.avg1(:,ch)=proj.sara.sim.avg1(:,ch)-mean(proj.sara.sim.avg1(:,ch));
    proj.sara.sim.avg1(:,ch)=proj.sara.sim.avg1(:,ch)/std(proj.sara.sim.avg1(:,ch));
    %find PSD
    [ob2.sara.ps2(:,ch),ob2.sara.freq2(:,ch)] = pwelch(proj.sara.sim.avg1(:,ch),fs,fs/2);
    subplot(8,2,ch)
    title(num2str(ch))
    plot(ob2.sara.freq2(:,ch)/(2*pi)*fs,ob2.sara.ps2(:,ch))
    xlim([0 40])
    % computing of the peak (crest) factor
    ob2.sara.rms2(ch) = sqrt(mean(ob2.sara.ps2(:,ch).^2));
    ob2.sara.peak2(ch) = max(abs(ob2.sara.ps2(:,ch)));
    ob2.sara.Q2(ch) = 20*log10(ob2.sara.peak2(ch)/ob2.sara.rms2(ch));
    disp(['Peak (crest) factor Q = ' num2str(ob2.sara.Q2(ch)) ' dB'])
    %spectrogram(ps1(:,ch), [], [], [], fs, 'yaxis')
end

%Identify a predictor for median nerve stimulation - upper limb
%each of the following have been identified as indicators for stimulation

%The earliest localized major component of the scalp recorded median 
%nerve SSEPs is the N20. It is recorded over the centroparietal region, 
%contralateral to the stimulated nerve and it represents the first 
%cortical response to the afferent somatosensory volley.

%N9 - average 11.8ms - window +- .87 ms - did not lead to good seperation
%N13 - average 13.9ms - window +- .47 ms - did not lead to good seperation
%N20 - average 21.5ms - window +- 1.20 ms - pretty good seperation
%N100 or N1 - average 100 - window +-20 ms - pretty good seperation
%
%Extract epochs from 80 and 120 milliseconds
times = 0:1/fs:1;
times = times(1:1200);
start = find(times >= .08,1);
finish = find(times >= .120,1);


for ch = 1:channels
[ERP.tc1v(ch),ERP.tc1i(ch)] = min(proj.tanya.consec.avg1(start:finish,ch));
[ERP.ts1v(ch),ERP.ts1i(ch)] = min(proj.tanya.sim.avg1(start:finish,ch));
[ERP.tc2v(ch),ERP.tc2i(ch)] = min(proj.tanya.consec.avg2(start:finish,ch));
[ERP.ts2v(ch),ERP.ts2i(ch)] = min(proj.tanya.sim.avg2(start:finish,ch));
[ERP.tc4v(ch),ERP.tc4i(ch)] = min(proj.tanya.consec.avg4(start:finish,ch));
[ERP.ts4v(ch),ERP.ts4i(ch)] = min(proj.tanya.sim.avg4(start:finish,ch));

[ERP.rc1v(ch),ERP.rc1i(ch)] = min(proj.roohi.consec.avg1(start:finish,ch));
[ERP.rs1v(ch),ERP.rs1i(ch)] = min(proj.roohi.sim.avg1(start:finish,ch));
[ERP.rc2v(ch),ERP.rc2i(ch)] = min(proj.roohi.consec.avg2(start:finish,ch));
[ERP.rs2v(ch),ERP.rs2i(ch)] = min(proj.roohi.sim.avg2(start:finish,ch));
[ERP.rc4v(ch),ERP.rc4i(ch)] = min(proj.roohi.consec.avg4(start:finish,ch));
[ERP.rs4v(ch),ERP.rs4i(ch)] = min(proj.roohi.sim.avg4(start:finish,ch));

[ERP.sc1v(ch),ERP.sc1i(ch)] = min(proj.sara.consec.avg1(start:finish,ch));
[ERP.ss1v(ch),ERP.ss1i(ch)] = min(proj.sara.sim.avg1(start:finish,ch));
[ERP.sc2v(ch),ERP.sc2i(ch)] = min(proj.sara.consec.avg2(start:finish,ch));
[ERP.ss2v(ch),ERP.ss2i(ch)] = min(proj.sara.sim.avg2(start:finish,ch));
[ERP.sc4v(ch),ERP.sc4i(ch)] = min(proj.sara.consec.avg4(start:finish,ch));
[ERP.ss4v(ch),ERP.ss4i(ch)] = min(proj.sara.sim.avg4(start:finish,ch));
end
predvc1 = [ERP.tc1v ERP.rc1v ERP.sc1v];
predvc2 = [ERP.tc2v ERP.rc2v ERP.sc2v];
predvc4 = [ERP.tc4v ERP.rc4v ERP.sc4v];

predvs1 = [ERP.ts1v ERP.rs1v ERP.ss1v];
predvs2 = [ERP.ts2v ERP.rs2v ERP.ss2v];
predvs4 = [ERP.ts4v ERP.rs4v ERP.ss4v];

predic1 = [ERP.tc1i ERP.rc1i ERP.sc1i];
predic2 = [ERP.tc2i ERP.rc2i ERP.sc2i];
predic4 = [ERP.tc4i ERP.rc4i ERP.sc4i];

predis1 = [ERP.ts1i ERP.rs1i ERP.ss1i];
predis2 = [ERP.ts2i ERP.rs2i ERP.ss2i];
predis4 = [ERP.ts4i ERP.rs4i ERP.ss4i];

state1 = [predvc1 ; predvs1];
state2 = [predvc2 ; predvs2];
state4 = [predvc4 ; predvs4];

%        1     2     3     4     5     6     7      8     9      10     11.  12     13     14    15     16
grp = {'Fpz'; 'Fz'; 'T7'; 'T8'; 'C3'; 'C4'; 'C5'; 'C6'; 'Cp3'; 'Cp4'; 'Cz'; 'Cpz'; 'Pz'; 'Po7'; 'Po8'; 'Oz'};

elec = [grp; grp; grp]';

stattes = [ones(16,1); ones(16,1)*2; 4*ones(16,1)];

R = table(state1',state2',state4',elec',stattes);
R.Properties.VariableNames = {'left','right','rest','channel','state'}


observations = [state1';state2';state4'];
channels = [elec';elec';elec'];
stattes = [ones(16*3,1); ones(16*3,1)*2; 4*ones(16*3,1)];

experiment = table(observations,channels,stattes);