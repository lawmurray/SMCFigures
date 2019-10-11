/**
 * Draw ancestry tree figures.
 */
program draw(width:Integer <- 1024, height:Integer <- 256, M:Integer <- 3,
    T:Integer <- 20, N:Integer <- 5, seed:Integer?) {
  if seed? {
    global.seed(seed!);
  }
    
  /* generate tree */
  W:Real[M,T,N];
  A:Integer[M,T,N];
  B:Integer[M,T];
  for auto m in 1..M {
    for auto n in 1..N {
      W[m,1,n] <- logpdf_gaussian(simulate_gaussian(0.0, 1.0), 0.0, 1.0);
    }
    A[m,1,1..N] <- iota(1, N);
    for auto t in 2..T {
      auto w <- W[m,t - 1,1..N];
      A[m,t,1..N] <- offspring_to_ancestors(simulate_multinomial(length(w), norm_exp(w)));
      for auto n in 1..N {
        W[m,t,n] <- logpdf_gaussian(simulate_gaussian(0.0, 1.0), 0.0, 1.0);
      }
      if m > 1 {
        A[m,t,B[m - 1,t]] <- B[m - 1,t - 1];
        W[m,t,B[m - 1,t]] <- W[m - 1,t,B[m - 1,t]];
      }
    }

    /* choose path */
    B[m,T] <- ancestor(W[m,T,1..N]);
    for auto r in 1..T-1 {
      auto t <- T - r;
      B[m,t] <- A[m,t + 1,B[m,t + 1]];
    }
  }

  /* draw frames */
  for auto m in 1..M {
    for auto t in 1..T {
      auto surface <- createPDF("figs/ancestry_" + m + "_" + t + ".pdf", width, height);
      auto cr <- create(surface);

      /* background */
      cr.setSourceRGB(0.95, 0.95, 0.95);
      cr.rectangle(0, 0, width, height);
      cr.fill();
  
      /* border and grid */
      cr.setSourceRGB(0.9, 0.9, 0.9);
      cr.setLineWidth(1);
      cr.rectangle(0, 0, width - 1, height - 1);
      for auto s in 1..T {
        cr.moveTo(x(s, T, width), 0);
        cr.lineTo(x(s, T, width), height);
      }
      cr.stroke();
  
      /* complete ancestry */
      cr.setSourceRGB(0.8, 0.8, 0.8);
      cr.setLineWidth(1);
      for auto s in 1..t {
        for auto n in 1..N {
          auto w <- exp(W[m,s,n]);
          auto a <- A[m,s,n];
        
          /* point for weighted particle */
          cr.arc(x(s, T, width), y(n, N, height), 10.0*sqrt(w), 0.0, 2.0*π);
          cr.fill();
          if s > 1 {
            /* line back to ancestor */
            cr.moveTo(x(s, T, width), y(n, N, height));
            cr.lineTo(x(s - 1, T, width), y(a, N, height));
            cr.stroke();
          }
        }
      }

      /* surviving ancestry */
      cr.setSourceRGB(0.0, 0.0, 0.0);
      cr.setLineWidth(1);
      for auto n in 1..N {
        auto a <- n;
        for auto r in 1..t {
          auto s <- t - r + 1;
          auto w <- exp(W[m,s,a]);
        
          /* point for weighted particle */
          cr.arc(x(s, T, width), y(a, N, height), 10.0*sqrt(w), 0.0, 2.0*π);
          cr.fill();
          if s > 1 {
            /* line back to ancestor */
            cr.moveTo(x(s, T, width), y(a, N, height));
            a <- A[m,s,a];
            cr.lineTo(x(s - 1, T, width), y(a, N, height));
            cr.stroke();
          }
        }
      }

      /* conditioned path */
      if m > 1 {
        cr.setSourceRGB(0.3373, 0.7059, 0.9137);
        cr.setLineWidth(3);
        for auto t in 2..T {
          auto b <- B[m - 1,t];
          auto a <- A[m - 1,t,b];
          auto w <- exp(W[m - 1,t,b]);
          
          cr.arc(x(t, T, width), y(b, N, height), 12.0*sqrt(w), 0.0, 2.0*π);
          cr.fill();
          cr.moveTo(x(t, T, width), y(b, N, height));
          cr.lineTo(x(t - 1, T, width), y(a, N, height));
          cr.stroke();
        }
      }
  
      /* destroy the surface (triggers save) */
      cr.destroy();
      surface.destroy();
    }
  }
}

function x(t:Integer, T:Integer, width:Integer) -> Integer {
  return Integer(width*(t - 0.5)/T);
}

function y(n:Integer, N:Integer, height:Integer) -> Integer {
  return Integer(height*(n - 0.5)/N);
}
