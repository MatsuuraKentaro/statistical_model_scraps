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
  vector[N] log_cumsum = log(reverse(cumulative_sum(reverse(exp(Xb)))));
  vector[N] pl = Event .* (Xb - log_cumsum);
  
  target += sum(pl);
  
  b[1:D] ~ std_normal();
}
