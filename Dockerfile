FROM gcc:4.9
WORKDIR /app
COPY . .
RUN gcc -o geneticAlgo geneticalgorithm.cpp -lstdc++
CMD ["./geneticAlgo"]
