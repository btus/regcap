# convert files from CRLF to LF
# Usage: rmcr filespec
#!/bin/bash
for f in "$@"
do
	cp "$f" /tmp/temp
        tr -d '\r' < /tmp/temp > "$f"
done
