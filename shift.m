function [t0, x0, u0] = shift(Td, t0, x0, u,f)

st = x0; con = u(1,:)';
f_value = f(st,con);

st = f_value;

x0 = full(st); % Novo inicijalno stanje sustava (iduci korak MPC petlje)

t0 = t0 + Td;

u0 = [u(2:size(u,1),:);u(size(u,1),:)]; % Inicijaliacija opt. vektora za iduci korak, koriste se izracunate vrijednosti prethodnog koraka
% Prva vrijednost (koristena) se makiva i za zadnje mjesto se ponavlja
% zadnja vrijednost prethodnog koraka

end