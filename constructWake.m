function [spWake,teWake] = constructWake(SP,TE,X,Y,strengths,Vinf,alpha,Nw)

spWake = generateStreamline(SP,true,2e-5,X,Y,strengths,Vinf,alpha,Nw);
teWake = generateStreamline(TE,true,2e-5,X,Y,strengths,Vinf,alpha,Nw);

end