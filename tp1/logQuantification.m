function [ newSound ] = logQuantification( sound, mu )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    newSound = zeros(length(sound),1)
    for x=1:length(sound)
        newSound(x) = sign(sound(x)) * ( log(1 + mu * sound(x)) / log(1 + mu));
    end
    
end

