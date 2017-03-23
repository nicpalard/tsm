s = audioread('samples/Toms_diner.wav');
s2 = s(1:2:length(s));
length(s2)
figure;
plot(s);
figure;
plot(s2);
%sound(s2, 22050);

s3 = resample(s2, 44100, 22050);
s4(:,1) = s; % Channel left
s4(:,2) = s3;% Channel right
sound(s4, 44100);

% Linear Quantification

soundQuant = logQuantification(s2, 255);
audiowrite(soundQuant, Fe, 22050, 8, 'Test.au');