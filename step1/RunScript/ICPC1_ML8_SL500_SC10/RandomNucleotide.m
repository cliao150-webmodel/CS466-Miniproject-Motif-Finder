function [Site,Count] = RandomNucleotide(rnum,profile)

Site = rnum;
Count = zeros(length(Site),4);

%   A
I1 = Site < sum(profile(1:1));
Site(I1) = 1;
Count(I1,1) = 1;

%   C
I2 = Site < sum(profile(1:2));
Site(I2) = 2;
Count(I2,2) = 1;

%   G
I3 = Site < sum(profile(1:3));
Site(I3) = 3;
Count(I3,3) = 1;

%   T
I4 = Site < 1.0;
Site(I4) = 4;
Count(I4,4) = 1;

end

