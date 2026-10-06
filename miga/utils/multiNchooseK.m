function val = multiNchooseK(a,g)
val = 1;
for i = 1:length(a)
    val = val * nchoosek(a(i),g(i));
end
end