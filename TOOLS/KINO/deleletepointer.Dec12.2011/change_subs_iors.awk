BEGIN{
 str=strftime("%b.%d.%Y");
 comment="ckino " str ": ";
}
/^ .*rv_p_orhoca/{
 print comment,$0
 if (match($0,"sspec\\(is\\)%rv_p_orhoc => rv_p_orhoca")){next}
 if (match($0,"rv_p_orhoca => sspec\\(is\\)%rv_p_orhoc")){next}
 if (match($0,"pointer *::")){next}
 gsub("rv_p_orhoca","sspec(is)%rv_p_orhoc")
}
/^ .*rv_p_oves/{
 print comment,$0
 if (match($0,"rv_p_oves => spot%rv_p_oves")){next}
 if (match($0,"spot%rv_p_oves => rv_p_oves")){next}
 if (match($0,"pointer *::")){next}
 gsub("rv_p_oves","spot%rv_p_oves")
}
/^ .*rv_p_ov1/{
 print comment,$0
 if (match($0,"ssite\\(ib\\)%rv_p_ov1 => rv_p_ov1")){next}
 if (match($0,"rv_p_ov1 => ssite\\(ib\\)%rv_p_ov1")){next}
 if (match($0,"pointer *::")){next}
 gsub("rv_p_ov1","ssite(ib)%rv_p_ov1")
}
/^ .*rv_p_ov0/{
 print comment,$0
 if (match($0,"ssite\\(ib\\)%rv_p_ov0 => rv_p_ov0")){next}
 if (match($0,"rv_p_ov0 => ssite\\(ib\\)%rv_p_ov0")){next}
 if (match($0,"pointer *::")){next}
 gsub("rv_p_ov0","ssite(ib)%rv_p_ov0")
}

{print}
