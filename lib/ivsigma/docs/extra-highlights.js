document.addEventListener('DOMContentLoaded', (event) => {
    for (let el of document.querySelectorAll('code.highlight.language-cpp')) {
        hljs.highlightBlock(el)
        el.style.display = "inline"
    }
});
