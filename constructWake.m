function [spWake,teWake] = constructWake(SP,TE,X,Y,strengths,Vinf,alpha,Nw)

spWake = generateStreamline(SP,X,Y,strengths,Vinf,alpha,Nw);
teWake = generateStreamline(TE,X,Y,strengths,Vinf,alpha,Nw);

end