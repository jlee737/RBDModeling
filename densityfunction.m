function density = densityfunction(dist, library_size)
    x = linspace(0, 300, 100000);
    pdfValues = pdf(dist, x);
    density = pdfValues * library_size;
end
