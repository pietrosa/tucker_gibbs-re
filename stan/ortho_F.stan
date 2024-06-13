data {
  int<lower=1> n; // dim of orthogonal vecs
  int<lower=1, upper=n> r; // # of orthogonal vecs
  matrix[n,r] F; // MvMVF parameter
}
parameters {
  matrix[n, r] X;
}
transformed parameters {
  matrix[n, r] Q = qr_thin_Q(X);
}
model {
  target += trace(-0.5*X'*X + F'*Q);
}
generated quantities {
  //
}
