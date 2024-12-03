Tx = 128;
Rx = 128;


H = randn(Rx,Tx);

x = randn(Tx,1);

n = randn(Rx,1);

[Q,R] = qr(H);

rx = R * x;

nz = Q' *n;


Es = sum(rx.^2);

En = sum(nz.^2);

snr = Es/En;