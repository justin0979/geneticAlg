FROM gcc:13
WORKDIR /app
COPY . .
RUN gcc -o main ./src/main.cpp -lstdc++
CMD ["./main"]
