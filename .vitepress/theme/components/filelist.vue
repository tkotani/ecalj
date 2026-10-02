<template>
  <ul>
    <li v-for="file in pdfFiles" :key="file.link">
      <a :href="file.link" target="_blank">{{ file.name }}</a>
    </li>
  </ul>
</template>

<script setup>
  const pdfFiles = Object.entries(import.meta.glob('/unchecked/**/*.pdf')).map(([path]) => {
  const name = path.split('/').pop(); // ファイル名を取得
  const link = `/ecaljdoc${path}`; // ecaljdoc をリンクに追加
  console.log(link); // 出力: /ecaljdoc/unchecked/file.txt
  return { link: link, name };
});
</script>
