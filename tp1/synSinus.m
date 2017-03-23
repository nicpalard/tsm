function [ tab ] = synSinus(size, f, amp, phi, Fech)
%SYN_SINUS Summary of this function goes here
%   Detailed explanation goes here
    i = 0:size;
    tab = amp * cos(2* pi * f * i/Fech + phi);
end
