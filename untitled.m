p = diag([2 1 2]);
a = [0 0 4; 2 0 0; 0 1 0];
l = 2*eye(3) - a;
q = p*l + l'*p;

% Check if PSD
eig(q)  % Should be [0, positive, positive]

% The key: negative off-diagonals don't prevent PSD!