FROM debian:bookworm-slim AS builder

RUN apt-get update && apt-get install -y \
    build-essential \
    git \
    && rm -rf /var/lib/apt/lists/*

COPY . /app

RUN cd /app && make all
FROM debian:bookworm-slim
COPY --from=builder /app/bin/flex_demux_mtx /app/bin/call_features /usr/local/bin/