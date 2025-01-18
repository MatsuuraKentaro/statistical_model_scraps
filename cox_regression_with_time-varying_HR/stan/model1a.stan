data {
  int N;
  int D;
  matrix[N, D] X;
  vector[N] Event;
}

parameters {
  vector[D] b;
}

model {
  vector[N] Xb = X * b;
  vector[N] pl;
  
  for (n in 1:N) {
    pl[n] = Event[n] * (Xb[n] - log_sum_exp(Xb[n:N]));
  }
  
  target += sum(pl);
  
  b[1:D] ~ std_normal();
}
