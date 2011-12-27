/^ .*rv_p_oqp/{
print FNR,$0
}
/^ .*if /{
print FNR,$0
}
/^ .*endif/{
print FNR,$0
}
/^ .*do /{
print FNR,$0
}
/^ .*enddo/{
print FNR,$0
}
/^ .*goto/{
print FNR,$0
}
/^ .[0-9]+/{
if (!match($0,"format")) {
print FNR,$0
}
}
