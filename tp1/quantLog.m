
N=8;
mu=255;
N_valeurs=2^(N-1);
s=wavread('Toms_diner.wav');

% compression
s_comp=sign(s).*(log(1+abs(s).*mu)/log(1+mu));

% quantification lineaire
s_quant=round(s_comp.*N_valeurs)./N_valeurs;

% expansion
s_exp=(sign(s_quant).*(1/mu)).*(exp(abs(s_quant).*(log(1+mu)))-1);

auwrite(s_exp,44100,8,'mu','tomsmu.au');
s8=wavread('Toms_diner_8bits.wav');
auwrite(s8,44100,8,'linear','tomslin.au');







